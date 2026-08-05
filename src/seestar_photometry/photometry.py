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


def _detect(plane, thresh):
    """Background-subtract a plane and detect sources. Returns (bkg, sub, objs)."""
    bkg = sep.Background(plane)
    data_sub = plane - bkg
    objs = sep.extract(data_sub, thresh, err=bkg.globalrms)
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


def _bright_round(data_sub, objs, rms, snr_min, ref_radius=10.0):
    """Mask of bright, round, unflagged sources for PSF/COG measurement."""
    flux, fluxerr, _ = sep.sum_circle(
        data_sub, objs["x"], objs["y"], ref_radius, err=rms
    )
    return (flux / fluxerr > snr_min) & (objs["b"] / objs["a"] > 0.7)


def _fwhm_from(data_sub, objs, rms, snr_min=SNR_PSF):
    """Median FWHM (px) from the half-light radii of bright round stars."""
    flux, _, _ = sep.sum_circle(data_sub, objs["x"], objs["y"], 10.0, err=rms)
    rhalf, flag = sep.flux_radius(
        data_sub, objs["x"], objs["y"], 6.0 * objs["a"], 0.5,
        normflux=flux, subpix=5,
    )
    sel = _bright_round(data_sub, objs, rms, snr_min) & (flag == 0)
    return GAUSS_FWHM * float(np.median(rhalf[sel]))


def _cog_from(data_sub, objs, rms, radii, snr_min=SNR_PSF):
    """Median curve of growth from bright, round, isolated stars.

    Isolation matters more than it looks: a neighbour inside the outermost aperture
    inflates the "total" flux and flattens the curve, which would push the measured
    90% radius outward and quietly enlarge the aperture for the whole frame.
    """
    x, y = objs["x"], objs["y"]
    measured = [sep.sum_circle(data_sub, x, y, r, err=rms) for r in radii]
    flux = np.array([m[0] for m in measured])
    flag = np.array([m[2] for m in measured])
    nn_dist = cKDTree(np.column_stack([x, y])).query(
        np.column_stack([x, y]), k=2
    )[0][:, 1]
    clean = (
        _bright_round(data_sub, objs, rms, snr_min)
        & (flag.sum(axis=0) == 0)
        & (nn_dist > 2 * radii.max())
    )
    norm = flux[:, clean] / flux[-1, clean]
    t = Table()
    t["radius"] = radii
    t["flux_frac"] = np.median(norm, axis=1)
    t["flux_frac_std"] = np.std(norm, axis=1)
    t.meta["n_stars"] = int(clean.sum())
    return t


def measure_fwhm(frame, band="G", thresh=2.0, snr_min=SNR_PSF):
    """Median PSF FWHM (pixels) of a frame, from bright round stars."""
    i = BANDS.index(band)
    bkg, data_sub, objs = _detect(frame.data[i], thresh)
    return _fwhm_from(data_sub, objs, bkg.globalrms, snr_min)


def curve_of_growth(frame, band="G", radii=None, snr_min=SNR_PSF, thresh=2.0):
    """Median curve of growth from bright, isolated stars in one band.

    Measures enclosed flux in a series of circular apertures, normalised to the
    outermost (~total). Returns a table with ``radius``, ``flux_frac``,
    ``flux_frac_std`` and ``meta["n_stars"]``.
    """
    radii = COG_RADII if radii is None else np.asarray(radii, dtype=float)
    i = BANDS.index(band)
    bkg, data_sub, objs = _detect(frame.data[i], thresh)
    return _cog_from(data_sub, objs, bkg.globalrms, radii, snr_min)


def aperture_for_enclosed_flux(cog, enclosed=0.90):
    """Radius (px) enclosing ``enclosed`` fraction of the flux, from a COG."""
    return float(np.interp(enclosed, cog["flux_frac"], cog["radius"]))


def aperture_correction(frame, radius, band="G", cog=None):
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
        cog = curve_of_growth(frame, band=band)
    frac = np.interp(radius, cog["radius"], cog["flux_frac"])
    return 2.5 * np.log10(frac)


def _size_aperture(data_sub, objs, rms, fwhm, frame_path, band,
                   aperture=None, enclosed=0.90, n_fwhm=None):
    """Pick a circular-aperture radius (px) for one band. Returns ``(radius, cog)``.

    The default sizes each band to the radius enclosing ``enclosed`` of its flux
    from a curve of growth, so the aperture is chromatic (matching the Seestar PSF)
    and its aperture correction is identical across bands. An explicit ``aperture``
    (shared radius) or ``n_fwhm`` (multiple of the band FWHM) overrides that. Falls
    back to ``1.2 * fwhm`` (~the 90% aperture) when a band has too few isolated
    stars for a COG.
    """
    cog = None
    if aperture is not None:
        ap = float(aperture)
    elif n_fwhm is not None:
        ap = float(n_fwhm * fwhm)
    else:
        cog = _cog_from(data_sub, objs, rms, COG_RADII)
        ap = aperture_for_enclosed_flux(cog, enclosed)
        if not np.isfinite(ap):
            ap = 1.2 * fwhm
    if not np.isfinite(ap):
        raise ValueError(f"cannot size {band}-band aperture for {frame_path}")
    return ap, cog


def _per_band_setup(frame, thresh, aperture, enclosed, n_fwhm):
    """Run detection and aperture sizing on all three planes.

    Yields ``(i, band, bkg, data_sub, objs, fwhm, ap, cog)`` per band. Shared by
    :func:`extract_sources` and :func:`forced_photometry` so the two can never
    drift apart in how they size an aperture -- if they did, a forced measurement
    would not be on the same photometric system as the zero point fitted from the
    detections.
    """
    for i, band in enumerate(BANDS):
        bkg, data_sub, objs = _detect(frame.data[i], thresh)
        fwhm = _fwhm_from(data_sub, objs, bkg.globalrms)
        ap, cog = _size_aperture(
            data_sub, objs, bkg.globalrms, fwhm, frame.path, band,
            aperture=aperture, enclosed=enclosed, n_fwhm=n_fwhm,
        )
        yield i, band, bkg, data_sub, objs, fwhm, ap, cog


def extract_sources(
    frame, thresh=2.0, aperture=None, enclosed=0.90, n_fwhm=None, bkgann=None
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
        frame, thresh, aperture, enclosed, n_fwhm
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
    aperture=None, enclosed=0.90, n_fwhm=None, bkgann=None,
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
    thresh, aperture, enclosed, n_fwhm, bkgann
        As in :func:`extract_sources`.

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
        frame, thresh, aperture, enclosed, n_fwhm
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
