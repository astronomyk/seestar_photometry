"""Source detection and fixed-aperture / forced photometry.

This is the standard photometry path. It uses a circular aperture sized per frame
**and per band** from each plane's point-spread function (default: the radius
enclosing 90% of the flux).

The Seestar PSF is chromatic: the cheap optics focus best near green, so both the
red and blue planes are measurably broader than green (typically R ~1.1x, B up to
~1.3x the green FWHM, varying frame to frame with focus). A single radius shared
across bands would therefore enclose a different flux fraction in each band and
bias the B-R colour the Gaia calibration relies on -- and because the broadening
tracks focus, that bias varies frame to frame and cannot be absorbed into a
constant colour term. Sizing each band to the same enclosed *fraction* instead
makes the aperture correction identical across bands, so it cancels in colours and
is absorbed into each per-frame, per-band zero point.

Kron/AUTO photometry lives in :mod:`kron`; it is kept only as a cross-check
(star/galaxy separation, total-flux comparison), not the standard path, because its
per-source adaptive aperture introduces source-dependent systematics that do not
cancel as cleanly.

Every entry point takes an optional ``mask`` -- a boolean array, ``True`` where pixels
should be ignored, following SEP's convention. It exists because **extended emission
corrupts the aperture, not the FWHM**. A field containing a bright nebula (M27, say)
feeds a crowd of resolved knots into the bright/round sample that sizes the aperture;
their growth curves never converge, so the median curve of growth is depressed at every
radius and :func:`aperture_for_enclosed_flux` slides far out -- measured on a real M27
frame, a 5 px green aperture became 19 px while the FWHM moved by 0.01 px, because a
median of per-source scalars shrugs off a minority that a median of whole *curves* does
not. Nothing detects this: the COG star count stays high, so the too-few-stars fallback
never fires and the bad aperture is used silently. ``mask`` is how you rule the region
out. See ``docs/photometry-design.md``.

See ``docs/photometry-design.md`` for the measurements behind these choices.
"""

from dataclasses import dataclass, field

import numpy as np
import sep
from astropy.table import Table, vstack
from scipy.spatial import cKDTree

from .frames import BANDS

#: FWHM of a Gaussian PSF in units of its half-light (50%-flux) radius.
GAUSS_FWHM = 2.0

#: SNR floor for a detection to enter the ``sources`` table. Deliberately the cut
#: that defines the sample -- the detection threshold is kept well below it.
SNR_MIN = 5.0

#: SNR a source must clear to be used for PSF / curve-of-growth measurement. High,
#: because these stars set the aperture for every other source in the frame.
SNR_PSF = 50.0

#: Radii (px) the curve of growth is sampled at. The outermost is taken as "total".
COG_RADII = np.arange(1.0, 21.0)

#: A curve-of-growth candidate is rejected when its semi-major axis exceeds this
#: multiple of the median semi-major axis of the frame's bright, round sources.
#:
#: Without a size cut the isolation criterion selects *for* extended objects in a
#: crowded field, which is the opposite of what it is for: only something large has no
#: close neighbour. On a real M27 frame (5469 green detections) the nearest-neighbour
#: distance was a median 10.6 px, so of 578 bright round sources exactly **one** cleared
#: the 40 px isolation cut -- M27 itself, at a = 42 px with b/a = 0.82. Roundness does
#: not exclude a nebula; nothing did. The aperture was then sized from the nebula's
#: growth curve, giving 19.0 px instead of 5.3.
COG_MAX_SIZE_RATIO = 3.0

#: Below this many stars a curve of growth is not trusted, and the aperture falls back
#: to ``1.2 * FWHM``. A COG built from one or two objects yields a *finite* radius, so
#: without an explicit floor the non-finite check alone lets it through silently.
MIN_COG_STARS = 5


@dataclass
class Extraction:
    """Result of a photometry run on one frame.

    Attributes
    ----------
    background : np.ndarray
        The ``(3, ny, nx)`` per-band background image.
    rms : np.ndarray
        Per-band global background RMS, shape ``(3,)`` in R, G, B order.
    sources : astropy.table.Table
        Photometry of all SNR > 5 detections, with a ``band`` column.
    aperture : np.ndarray or None
        Per-band circular aperture radius used (px), shape ``(3,)``.
    fwhm : np.ndarray or None
        Per-band median FWHM measured for the frame (px), shape ``(3,)``.
    cogs : list or None
        Per-band curve-of-growth tables, in R, G, B order (``None`` for a band
        sized without one). Kept for the diagnostic figures.
    frame : SeestarFrame or None
        The frame this came from (used to locate the WCS cache and read solver
        hints).
    """

    background: np.ndarray
    rms: np.ndarray
    sources: Table
    aperture: np.ndarray = None
    fwhm: np.ndarray = None
    cogs: list = field(default=None, repr=False)
    frame: object = field(default=None, repr=False)

    def solve_wcs(self, api_key=None, force=False, solver="nova", catalogue=None):
        """Solve (or load cached) the frame's WCS from this source catalogue.

        Reuses the green-band detections as the solver's source list, so no second
        extraction is needed. The solution is cached as a ``.wcs`` sidecar next to
        the FITS file. Returns an :class:`astropy.wcs.WCS`.

        ``catalogue`` is the reference table ``solver="local"`` pairs against, and is
        ignored by the other solvers.
        """
        from . import astrometry

        if self.frame is None:
            raise ValueError("Extraction has no frame; cannot solve/cache WCS")
        g = self.sources[self.sources["band"] == "G"]
        order = np.argsort(np.asarray(g["flux"]))[::-1]  # brightest first
        return astrometry.solve_from_sources(
            self.frame, np.asarray(g["x"])[order], np.asarray(g["y"])[order],
            api_key=api_key, force=force, solver=solver, catalogue=catalogue,
        )

    def match_gaia(self, gaia, wcs=None, tol_arcsec=2.0):
        """Cross-match this extraction's sources to a reference catalogue.

        Solves/loads the WCS (or uses the one passed), subsets ``gaia`` to this
        frame's footprint, and augments ``.sources`` **in place** with ``ra``,
        ``dec``, ``sep_arcsec`` and the matched catalogue columns (synthetic V,
        B/R, ``phot_variable_flag``, ...), masked where there is no match. Returns
        the augmented table as well.
        """
        from . import catalogs

        if self.frame is None:
            raise ValueError("Extraction has no frame; cannot cross-match")
        if wcs is None:
            wcs = self.solve_wcs()
        in_frame = catalogs.sources_in_frame(gaia, wcs, self.frame)
        self.sources = catalogs.crossmatch_table(self.sources, wcs, in_frame, tol_arcsec)
        return self.sources

    def band(self, band):
        """The ``sources`` rows for one band."""
        return self.sources[np.asarray(self.sources["band"]) == band]


def _as_mask(mask, shape):
    """Validate a mask and return it as a contiguous bool array, or ``None``."""
    if mask is None:
        return None
    mask = np.ascontiguousarray(mask, dtype=bool)
    if mask.shape != tuple(shape):
        raise ValueError(f"mask shape {mask.shape} does not match plane {tuple(shape)}")
    return mask


def _centroid_masked(objs, mask):
    """Which detections have their centroid on a masked pixel."""
    ny, nx = mask.shape
    x = np.clip(np.rint(np.asarray(objs["x"], dtype=float)), 0, nx - 1)
    y = np.clip(np.rint(np.asarray(objs["y"], dtype=float)), 0, ny - 1)
    bad = ~(np.isfinite(objs["x"]) & np.isfinite(objs["y"]))
    return mask[y.astype(int), x.astype(int)] | bad


def sky_mask(shape, wcs, ra, dec, radius_arcsec):
    """Boolean mask (``True`` = ignore) of circular sky regions, for ``mask=``.

    The natural way to specify an exclusion is on the sky -- "M27 is at this RA/Dec and
    about this big" -- not in pixels, which move frame to frame under dithering and
    Alt-Az field rotation. Radii are converted with the *solved* plate scale.

    Parameters
    ----------
    shape : tuple
        ``(ny, nx)`` of the plane being masked. Take it from ``frame.shape``, never a
        header.
    wcs : astropy.wcs.WCS
    ra, dec : float or array-like
        Region centres, degrees.
    radius_arcsec : float or array-like
        Region radii, arcseconds. Broadcast against ``ra``/``dec``.

    Returns
    -------
    np.ndarray
        Boolean ``(ny, nx)``, ``True`` inside any region.
    """
    from astropy.wcs.utils import proj_plane_pixel_scales

    ny, nx = shape
    ra = np.atleast_1d(np.asarray(ra, dtype=float))
    dec = np.atleast_1d(np.asarray(dec, dtype=float))
    radius = np.broadcast_to(np.atleast_1d(np.asarray(radius_arcsec, dtype=float)),
                            ra.shape)
    scale = float(np.mean(proj_plane_pixel_scales(wcs.celestial))) * 3600.0
    x, y = wcs.world_to_pixel_values(ra, dec)
    yy, xx = np.indices((ny, nx))
    mask = np.zeros((ny, nx), dtype=bool)
    for cx, cy, r in zip(np.atleast_1d(x), np.atleast_1d(y), radius):
        if not (np.isfinite(cx) and np.isfinite(cy)):
            continue
        r_px = float(r) / scale
        mask |= (xx - float(cx)) ** 2 + (yy - float(cy)) ** 2 < r_px ** 2
    return mask


def _detect(plane, thresh, mask=None):
    """Background-subtract a plane and detect sources. Returns (bkg, sub, objs).

    With a ``mask``, the background is estimated ignoring masked pixels and detections
    whose **centroid** lands on one are dropped.

    The mask is deliberately *not* handed to ``sep.extract``. Doing so does not remove a
    bright extended source: SEP excludes the masked pixels from the footprint but the
    surrounding flux still clears the threshold, so one nebula fragments into a ring of
    detections around the mask boundary -- whose centroids are all *outside* the mask and
    would survive any centroid cut. Rejecting on the centroid of the un-masked
    segmentation removes the object as one piece instead.
    """
    mask = _as_mask(mask, plane.shape)
    bkg = sep.Background(plane, mask=mask)
    data_sub = plane - bkg
    objs = sep.extract(data_sub, thresh, err=bkg.globalrms)
    if mask is not None and len(objs):
        objs = objs[~_centroid_masked(objs, mask)]
    return bkg, data_sub, objs


def _max_in_aperture(data, x, y, r):
    """Highest raw pixel value within radius ``r`` of each ``(x, y)``.

    Computed on the raw (non background-subtracted) plane so it can be compared
    directly against the 16-bit ceiling for saturation masking.
    """
    ny, nx = data.shape
    rr = int(np.ceil(r))
    yy, xx = np.ogrid[-rr:rr + 1, -rr:rr + 1]
    circle = xx ** 2 + yy ** 2 <= r ** 2
    out = np.empty(len(x))
    for k in range(len(x)):
        if not (np.isfinite(x[k]) and np.isfinite(y[k])):
            out[k] = np.nan  # degenerate detection with a NaN centroid
            continue
        cx, cy = int(round(x[k])), int(round(y[k]))
        y0, y1 = max(cy - rr, 0), min(cy + rr + 1, ny)
        x0, x1 = max(cx - rr, 0), min(cx + rr + 1, nx)
        cmask = circle[y0 - (cy - rr):y1 - (cy - rr), x0 - (cx - rr):x1 - (cx - rr)]
        region = data[y0:y1, x0:x1][cmask]
        out[k] = region.max() if region.size else np.nan
    return out


def _bright_round(data_sub, objs, rms, snr_min, ref_radius=10.0, mask=None):
    """Mask of bright, round, unflagged sources for PSF/COG measurement.

    A ``mask`` additionally rejects sources whose reference aperture *overlaps* a masked
    pixel, via SEP's ``APER_HASMASKED`` flag -- a knot just outside the region is as bad
    for the aperture as one inside it. The bit cannot be set without a mask, so passing
    ``None`` leaves the selection bit-identical to before.
    """
    flux, fluxerr, flag = sep.sum_circle(
        data_sub, objs["x"], objs["y"], ref_radius, err=rms, mask=mask
    )
    return (
        (flux / fluxerr > snr_min)
        & (objs["b"] / objs["a"] > 0.7)
        & ((flag & sep.APER_HASMASKED) == 0)
    )


def _fwhm_from(data_sub, objs, rms, snr_min=SNR_PSF, mask=None):
    """Median FWHM (px) from the half-light radii of bright round stars."""
    flux, _, _ = sep.sum_circle(data_sub, objs["x"], objs["y"], 10.0, err=rms, mask=mask)
    rhalf, flag = sep.flux_radius(
        data_sub, objs["x"], objs["y"], 6.0 * objs["a"], 0.5,
        normflux=flux, subpix=5,
    )
    sel = _bright_round(data_sub, objs, rms, snr_min, mask=mask) & (flag == 0)
    return GAUSS_FWHM * float(np.median(rhalf[sel]))


def _cog_from(data_sub, objs, rms, radii, snr_min=SNR_PSF, mask=None,
              size_ratio=COG_MAX_SIZE_RATIO, isolation=None):
    """Median curve of growth from bright, round, isolated, **point-like** stars.

    Isolation matters more than it looks: a neighbour inside the outermost aperture
    inflates the "total" flux and flattens the curve, which would push the measured
    90% radius outward and quietly enlarge the aperture for the whole frame.

    But isolation alone selects *for* extended objects in a crowded field -- only
    something large has no close neighbour -- and roundness does not exclude a nebula.
    Hence the size cut (``size_ratio``, see :data:`COG_MAX_SIZE_RATIO`), which is what
    actually keeps a planetary nebula or galaxy out of the sample. ``mask`` handles the
    complementary case of emission you know about a priori: masked pixels anywhere in a
    star's apertures set ``APER_HASMASKED``, and the all-radii flag check drops it.

    ``isolation`` overrides the nearest-neighbour distance (px) a candidate must clear,
    default ``2 * radii.max()``. In a rich field that default can admit almost nothing --
    on a real M27 frame the median nearest-neighbour distance was 10.6 px against a 40 px
    cut -- so lowering it is how you get a measured COG rather than the FWHM fallback.
    Lower it knowingly: the flux at the outer radii is what normalises the curve.
    """
    x, y = objs["x"], objs["y"]
    if not len(objs):
        return _empty_cog(radii)
    measured = [sep.sum_circle(data_sub, x, y, r, err=rms, mask=mask) for r in radii]
    flux = np.array([m[0] for m in measured])
    flag = np.array([m[2] for m in measured])
    nn_dist = cKDTree(np.column_stack([x, y])).query(
        np.column_stack([x, y]), k=2
    )[0][:, 1]
    bright = _bright_round(data_sub, objs, rms, snr_min, mask=mask)
    isolation = 2 * radii.max() if isolation is None else float(isolation)
    clean = (
        bright
        & _point_like(objs, bright, size_ratio)
        & (flag.sum(axis=0) == 0)
        & (nn_dist > isolation)
    )
    if not clean.any():
        return _empty_cog(radii)
    norm = flux[:, clean] / flux[-1, clean]
    t = Table()
    t["radius"] = radii
    t["flux_frac"] = np.median(norm, axis=1)
    t["flux_frac_std"] = np.std(norm, axis=1)
    t.meta["n_stars"] = int(clean.sum())
    return t


def _empty_cog(radii):
    """A curve of growth with no usable stars: all-nan, ``n_stars = 0``.

    Returned explicitly rather than letting ``np.median`` of an empty slice produce the
    same nans plus a RuntimeWarning, so an empty sample is a stated outcome the caller
    can branch on rather than a warning to be filtered out of the logs.
    """
    t = Table()
    t["radius"] = np.asarray(radii, dtype=float)
    t["flux_frac"] = np.full(len(radii), np.nan)
    t["flux_frac_std"] = np.full(len(radii), np.nan)
    t.meta["n_stars"] = 0
    return t


def _point_like(objs, reference, size_ratio=COG_MAX_SIZE_RATIO):
    """Which sources are no larger than ``size_ratio`` x the typical bright source.

    The scale comes from the median semi-major axis of ``reference`` (the frame's bright,
    round sources) rather than from the FWHM, so it needs no second pass and adapts to
    focus. A galaxy or planetary nebula sits orders of magnitude above it; a star sits at
    1 by construction.
    """
    a = np.asarray(objs["a"], dtype=float)
    if reference is None or not np.any(reference):
        return np.ones(len(a), dtype=bool)
    a_ref = float(np.median(a[reference]))
    if not np.isfinite(a_ref) or a_ref <= 0:
        return np.ones(len(a), dtype=bool)
    return a <= size_ratio * a_ref


def measure_fwhm(frame, band="G", thresh=2.0, snr_min=SNR_PSF, mask=None):
    """Median PSF FWHM (pixels) of a frame, from bright round stars.

    ``mask`` (``True`` = ignore) excludes a region such as a nebula. The FWHM is far more
    robust to extended emission than the aperture is -- it is a median over hundreds of
    per-source half-light radii, and on a real M27 frame masking moved it by 0.01 px -- so
    passing a mask here is a refinement, not a fix. :func:`curve_of_growth` is the one
    that needs it.
    """
    i = BANDS.index(band)
    bkg, data_sub, objs = _detect(frame.data[i], thresh, mask=mask)
    return _fwhm_from(data_sub, objs, bkg.globalrms, snr_min, mask=mask)


def curve_of_growth(frame, band="G", radii=None, snr_min=SNR_PSF, thresh=2.0,
                    mask=None, size_ratio=COG_MAX_SIZE_RATIO, isolation=None):
    """Median curve of growth from bright, isolated, point-like stars in one band.

    Measures enclosed flux in a series of circular apertures, normalised to the
    outermost (~total). Returns a table with ``radius``, ``flux_frac``,
    ``flux_frac_std`` and ``meta["n_stars"]``.

    **Check ``meta["n_stars"]``.** A curve built from one or two objects still returns a
    finite radius, and in a crowded or nebulous field that is the likely outcome -- see
    :data:`COG_MAX_SIZE_RATIO` and :data:`MIN_COG_STARS`. ``flux_frac`` is all-nan with
    ``n_stars = 0`` when nothing qualifies.

    ``mask`` (``True`` = ignore) rules out a region of extended emission; ``isolation``
    relaxes the nearest-neighbour requirement for crowded fields.
    """
    radii = COG_RADII if radii is None else np.asarray(radii, dtype=float)
    i = BANDS.index(band)
    bkg, data_sub, objs = _detect(frame.data[i], thresh, mask=mask)
    return _cog_from(data_sub, objs, bkg.globalrms, radii, snr_min, mask=mask,
                     size_ratio=size_ratio, isolation=isolation)


def aperture_for_enclosed_flux(cog, enclosed=0.90):
    """Radius (px) enclosing ``enclosed`` fraction of the flux, from a COG."""
    return float(np.interp(enclosed, cog["flux_frac"], cog["radius"]))


def aperture_correction(frame, radius, band="G", cog=None, mask=None):
    """Magnitude offset from a fixed-aperture mag to a total mag.

    Reads the enclosed-flux fraction at ``radius`` off the curve of growth and
    returns ``2.5 * log10(frac)`` (negative; the total flux is brighter). Add it to
    a circular-aperture instrumental magnitude to recover the total magnitude. Pass
    a precomputed ``cog`` to avoid rebuilding it.

    Note that the standard path never needs this: sizing every band to the same
    enclosed fraction makes the correction a constant that is absorbed into the
    zero point and cancels in colours. It is here for putting a total magnitude on
    an absolute scale.
    """
    if cog is None:
        cog = curve_of_growth(frame, band=band, mask=mask)
    frac = np.interp(radius, cog["radius"], cog["flux_frac"])
    return 2.5 * np.log10(frac)


def _size_aperture(data_sub, objs, rms, fwhm, frame_path, band,
                   aperture=None, enclosed=0.90, n_fwhm=None, mask=None,
                   isolation=None):
    """Pick a circular-aperture radius (px) for one band. Returns ``(radius, cog)``.

    The default sizes each band to the radius enclosing ``enclosed`` of its flux
    from a curve of growth, so the aperture is chromatic (matching the Seestar PSF)
    and its aperture correction is identical across bands. An explicit ``aperture``
    (shared radius) or ``n_fwhm`` (multiple of the band FWHM) overrides that.

    Falls back to ``1.2 * fwhm`` (~the 90% aperture) when the COG is not trustworthy:
    fewer than :data:`MIN_COG_STARS` stars, or a non-finite radius. The star-count floor
    is the load-bearing half -- a COG from a single object returns a *finite* radius, so
    the non-finite check alone would pass 19 px through as if it were measured.
    """
    cog = None
    if aperture is not None:
        ap = float(aperture)
    elif n_fwhm is not None:
        ap = float(n_fwhm * fwhm)
    else:
        cog = _cog_from(data_sub, objs, rms, COG_RADII, mask=mask,
                        isolation=isolation)
        ap = aperture_for_enclosed_flux(cog, enclosed)
        if cog.meta["n_stars"] < MIN_COG_STARS or not np.isfinite(ap):
            ap = 1.2 * fwhm
    if not np.isfinite(ap):
        raise ValueError(f"cannot size {band}-band aperture for {frame_path}")
    return ap, cog


def _per_band_setup(frame, thresh, aperture, enclosed, n_fwhm, mask=None,
                    isolation=None):
    """Run detection and aperture sizing on all three planes.

    Yields ``(i, band, bkg, data_sub, objs, fwhm, ap, cog)`` per band. Shared by
    :func:`extract_sources` and :func:`forced_photometry` so the two can never
    drift apart in how they size an aperture -- if they did, a forced measurement
    would not be on the same photometric system as the zero point fitted from the
    detections.
    """
    for i, band in enumerate(BANDS):
        bkg, data_sub, objs = _detect(frame.data[i], thresh, mask=mask)
        fwhm = _fwhm_from(data_sub, objs, bkg.globalrms, mask=mask)
        ap, cog = _size_aperture(
            data_sub, objs, bkg.globalrms, fwhm, frame.path, band,
            aperture=aperture, enclosed=enclosed, n_fwhm=n_fwhm, mask=mask,
            isolation=isolation,
        )
        yield i, band, bkg, data_sub, objs, fwhm, ap, cog


def extract_sources(
    frame, thresh=2.0, aperture=None, enclosed=0.90, n_fwhm=None, bkgann=None,
    mask=None, isolation=None,
):
    """Detect sources and measure fixed circular-aperture flux across R, G, B.

    The aperture is sized **per frame and per band** from each plane's PSF, so every
    band encloses the same flux fraction and its aperture correction is identical --
    keeping the B-R colour unbiased despite the chromatic Seestar PSF (see the
    module docstring). Only sources with SNR > 5 are kept.

    Parameters
    ----------
    thresh : float
        Detection threshold in units of the background RMS. Kept low so the
        SNR > 5 cut, not the detection step, defines the final sample.
    aperture : float, optional
        Explicit radius in pixels, shared by all bands. Overrides the per-band
        sizing -- use only if you specifically want a common radius.
    enclosed : float
        If sizing by curve of growth (the default), each band uses the radius
        enclosing this fraction of its flux (~1.2 x FWHM for ``0.90``).
    n_fwhm : float, optional
        If given (and ``aperture`` is not), size each band's aperture as
        ``n_fwhm * FWHM`` of that band instead of by enclosed flux.
    bkgann : tuple of float, optional
        ``(r_in, r_out)`` local-background annulus in pixels; if given, a
        per-source local background is subtracted.
    mask : np.ndarray, optional
        Boolean ``(ny, nx)``, ``True`` where pixels should be ignored (SEP's
        convention); build one with :func:`sky_mask`. Excludes the region from the
        background estimate, from detection, and from the aperture-sizing sample. It does
        **not** alter the flux of any source that is measured: a star whose aperture
        overlaps the mask is still summed over all its pixels, because silently changing a
        reported flux is worse than reporting a contaminated one.
    isolation : float, optional
        Nearest-neighbour distance (px) a curve-of-growth candidate must clear. Default
        ``2 * COG_RADII.max()``; lower it in crowded fields, where the default can leave
        too few stars and force the FWHM fallback.

    Returns
    -------
    Extraction
    """
    background = np.empty_like(frame.data)
    rms = np.empty(3, dtype=np.float32)
    apertures = np.empty(3)
    fwhms = np.empty(3)
    cogs = [None, None, None]
    tables = []
    for i, band, bkg, data_sub, objs, fwhm, ap, cog in _per_band_setup(
        frame, thresh, aperture, enclosed, n_fwhm, mask=mask, isolation=isolation
    ):
        background[i] = bkg.back()
        rms[i] = bkg.globalrms
        fwhms[i] = fwhm
        apertures[i] = ap
        cogs[i] = cog
        flux, fluxerr, flag = sep.sum_circle(
            data_sub, objs["x"], objs["y"], ap, err=bkg.globalrms, bkgann=bkgann,
        )
        snr = flux / fluxerr
        keep = snr > SNR_MIN
        max_pix = _max_in_aperture(frame.data[i], objs["x"], objs["y"], ap)
        t = Table()
        t["band"] = np.full(int(keep.sum()), band)
        t["x"] = objs["x"][keep]
        t["y"] = objs["y"][keep]
        t["flux"] = flux[keep]
        t["fluxerr"] = fluxerr[keep]
        t["snr"] = snr[keep]
        t["flag"] = flag[keep]
        t["max_pix_value"] = max_pix[keep]
        t["a"] = objs["a"][keep]
        t["b"] = objs["b"][keep]
        t["theta"] = objs["theta"][keep]
        tables.append(t)
    return Extraction(
        background=background, rms=rms, sources=vstack(tables),
        aperture=apertures, fwhm=fwhms, cogs=cogs, frame=frame,
    )


def forced_photometry(
    frame, ra, dec, wcs, source_id=None, thresh=2.0,
    aperture=None, enclosed=0.90, n_fwhm=None, bkgann=None, mask=None,
    isolation=None,
):
    """Forced fixed-aperture photometry at given sky positions, across R, G, B.

    Unlike :func:`extract_sources`, this does **not** detect: it measures the
    aperture flux at every ``(ra, dec)`` position (mapped to pixels via ``wcs``),
    whether or not a source was detected there. That gives a complete, consistent
    sample for light curves -- a star keeps a row in every frame, so series never go
    ragged and the comparison ensemble is identical frame to frame. The aperture is
    sized exactly as in :func:`extract_sources`, so forced fluxes sit on the same
    photometric system as the zero point.

    Parameters
    ----------
    frame : SeestarFrame
    ra, dec : array-like
        Sky positions (deg) to measure -- typically in-footprint catalogue sources.
    wcs : astropy.wcs.WCS
        The frame's WCS, to map ``(ra, dec)`` to pixels.
    source_id : array-like, optional
        A stable identifier per position (e.g. Gaia ``source_id``), carried through
        to the output table.
    thresh, aperture, enclosed, n_fwhm, bkgann, mask, isolation
        As in :func:`extract_sources`. Note that ``mask`` governs the background and the
        aperture sizing only: **every** requested position is still measured, including
        one that falls inside the mask, because forced photometry's contract is to return
        a row for every position asked for. If your target sits on masked emission, the
        flux is the aperture sum including that emission -- see :mod:`contamination` for
        estimating and removing a host/nebula component.

    Returns
    -------
    astropy.table.Table
        One row per ``(position, band)`` with ``source_id`` (if given), ``band``,
        ``x``, ``y``, ``flux``, ``fluxerr``, ``snr``, ``flag``, ``max_pix_value``
        and ``on_chip``. Per-band ``aperture`` and ``fwhm`` land in ``meta``.
    """
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    x, y = wcs.world_to_pixel_values(ra, dec)
    ny, nx = frame.shape
    on_chip = (x >= 0) & (x < nx) & (y >= 0) & (y < ny)

    apertures = np.empty(3)
    fwhms = np.empty(3)
    tables = []
    for i, band, bkg, data_sub, objs, fwhm, ap, _cog in _per_band_setup(
        frame, thresh, aperture, enclosed, n_fwhm, mask=mask, isolation=isolation
    ):
        fwhms[i] = fwhm
        apertures[i] = ap
        flux, fluxerr, flag = sep.sum_circle(
            data_sub, x, y, ap, err=bkg.globalrms, bkgann=bkgann,
        )
        max_pix = _max_in_aperture(frame.data[i], x, y, ap)
        t = Table()
        if source_id is not None:
            t["source_id"] = np.asarray(source_id)
        t["band"] = np.full(len(x), band)
        t["x"] = x
        t["y"] = y
        t["flux"] = flux
        t["fluxerr"] = fluxerr
        # An aperture entirely off the frame has zero error; nan is the right answer,
        # not a warning, since forced photometry measures every position by design.
        with np.errstate(invalid="ignore", divide="ignore"):
            t["snr"] = np.where(fluxerr > 0, flux / fluxerr, np.nan)
        t["flag"] = flag
        t["max_pix_value"] = max_pix
        t["on_chip"] = on_chip
        tables.append(t)
    out = vstack(tables)
    out.meta["aperture"] = apertures
    out.meta["fwhm"] = fwhms
    return out


def fit_background(extraction, band="G"):
    """Fit a 2nd-order 2D polynomial to a band's SEP background map.

    Returns a dict with the ``background`` map, the polynomial ``model`` and
    ``residual`` images, the ``pedestal`` (median background, ADU), the six
    ``coeffs`` (term order ``1, x, y, x^2, xy, y^2`` in normalised coords [-1, 1])
    and ``resid_std``. The polynomial captures any smooth gradient; localised
    bright-star halos remain in the residual.
    """
    bg = extraction.background[BANDS.index(band)].astype(float)
    ny, nx = bg.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    xn, yn = xx / nx * 2 - 1, yy / ny * 2 - 1
    design = np.column_stack([
        np.ones(bg.size), xn.ravel(), yn.ravel(),
        (xn ** 2).ravel(), (xn * yn).ravel(), (yn ** 2).ravel(),
    ])
    coeffs, *_ = np.linalg.lstsq(design, bg.ravel(), rcond=None)
    model = (design @ coeffs).reshape(ny, nx)
    residual = bg - model
    return {
        "background": bg,
        "model": model,
        "residual": residual,
        "pedestal": float(np.median(bg)),
        "coeffs": coeffs,
        "resid_std": float(residual.std()),
    }


def instrumental_mag(flux):
    """``-2.5 * log10(flux)``, nan where the flux is non-positive.

    Non-positive forced fluxes are normal, not an error: a faint or absent source
    measured in a fixed aperture on a background-subtracted image scatters either
    side of zero. They must become nan rather than raise or warn, because the
    light-curve tables measure every catalogue position in every frame.
    """
    flux = np.asarray(flux, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(flux > 0, -2.5 * np.log10(flux), np.nan)


def mag_error(snr):
    """Photometric error (mag) from SNR: ``1.0857 / snr``, nan for snr <= 0."""
    snr = np.asarray(snr, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(snr > 0, 1.0857 / snr, np.nan)
