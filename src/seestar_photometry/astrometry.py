"""Per-frame WCS, cached as a ``.wcs`` sidecar next to the FITS file.

The on-board Seestar WCS is **not** photometric-grade: positions are off by tens of
pixels (~1 arcmin), which drops the catalogue match rate to a few percent and makes
forced photometry meaningless. Every frame is therefore re-solved, and the result is
cached as a header-only FITS sidecar so the (expensive) solve happens exactly once
per frame and is shared between projects.

Several solvers, same sidecar, chosen with ``solver=``:

``"local"``
    Anchored on a reference catalogue rather than on an index of quads. **Not a blind
    solver, and it does not need to be**: the on-board WCS is wrong by ~1 arcmin, not
    by degrees, so the pointing is a perfectly good seed and what is left is a
    refinement. Needs no external binary, no network and no index files -- only the
    catalogue the pipeline already has -- so it is the one solver that works
    everywhere. See :func:`solve_local`.

``"astap"`` (default)
    A local, fully offline blind plate solver. Does its own star detection, handles
    the field rotation of Alt-Az frames, and takes ~a second per frame. Still the
    default because it needs no catalogue at all.

``"nova"``
    astrometry.net's web service, fed our own SEP source list (not the image) with
    the header pointing and pixel scale as hints. Useful as a fallback when ASTAP
    cannot solve a shallow stack, and it is what the historical results were built
    with. It intermittently drops connections under load, so callers should retry.

``lift`` covers the third case: a frame that already carries a *trustworthy*
solution in its header (CrowdSky plate-solves server-side and sets ``PLTSOLVD``), so
there is nothing to solve -- just copy it into a sidecar.

Degenerate (non-celestial) solutions are never cached, and a cached one that turns
out degenerate is transparently re-solved.

The astrometry.net API key is read from ``ASTROMETRY_KEY`` and never stored in the
repo. See ``docs/astrometry-and-gaia.md``.
"""

import os
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

#: Default location of the ASTAP command-line solver.
ASTAP_EXE = r"C:\Program Files\astap\astap_cli.exe"


def wcs_cache_path(frame_path):
    """Path of the ``.wcs`` sidecar for a given FITS path.

    A compression suffix is stripped first, so ``frame.fits.gz`` caches to
    ``frame.wcs`` rather than the surprising ``frame.fits.wcs`` that a naive suffix
    replacement gives. Gzipped frames are common -- astropy reads them transparently,
    so users do compress archives -- and a frame whose sidecar could not be found
    would be silently re-solved every run.
    """
    path = Path(frame_path)
    if path.suffix.lower() in (".gz", ".bz2", ".fz", ".zip"):
        path = path.with_suffix("")
    return path.with_suffix(".wcs")


def header_pixel_scale(header):
    """Nominal pixel scale (arcsec/pix) from the optics keywords."""
    return 206.265 * header["XPIXSZ"] / header["FOCALLEN"]


def pixel_scale(frame):
    """Pixel scale (arcsec/pix) for a frame, preferring a solved value.

    A plate-solved scale (``PIXSCALE``/``PLTSCALE``, set by CrowdSky) is measured
    rather than nominal, so it is used when present; otherwise it is derived from
    ``XPIXSZ``/``FOCALLEN``. The two agree to ~1% on the S50 (2.39 vs 2.37).
    """
    h = frame.header
    for key in ("PIXSCALE", "PLTSCALE"):
        value = h.get(key)
        if value is not None and np.isfinite(float(value)) and float(value) > 0:
            return float(value)
    return header_pixel_scale(h)


def _read_wcs(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return WCS(fits.getheader(path))


def _write_wcs(header, cache):
    """Cache a solution, refusing to store a degenerate one."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wcs = WCS(header)
    if not wcs.has_celestial:
        raise RuntimeError("refusing to cache a non-celestial WCS")
    fits.PrimaryHDU(header=header).writeto(cache, overwrite=True)
    return wcs


def load_wcs(frame):
    """The cached WCS for a frame, or ``None``.

    ``None`` means the frame is unsolved *or* the cached solution is degenerate
    (non-celestial), so a bad cache re-solves rather than poisoning the photometry.
    """
    path = wcs_cache_path(frame.path)
    if not path.exists():
        return None
    wcs = _read_wcs(path)
    return wcs if wcs.has_celestial else None


def has_wcs(frame_path):
    """Whether a usable sidecar exists, without loading the frame.

    Used by the batch runner to decide what still needs solving; taking a path
    rather than a frame avoids reading a few thousand FITS cubes just to skip them.
    """
    path = wcs_cache_path(frame_path)
    return path.exists() and _read_wcs(path).has_celestial


def lift(frame, force=False):
    """Copy a trustworthy header WCS into the sidecar, without solving.

    Only for frames whose header solution is known good -- CrowdSky sets
    ``PLTSOLVD = T`` after its own plate solve. Raises if the frame carries no such
    marker, because lifting the *Seestar's* on-board WCS would look like success and
    silently produce ~1 arcmin position errors.
    """
    cache = wcs_cache_path(frame.path)
    if not force and cache.exists() and _read_wcs(cache).has_celestial:
        return _read_wcs(cache)
    if not frame.header.get("PLTSOLVD"):
        raise RuntimeError(
            f"{frame.path} has no trustworthy header WCS (PLTSOLVD not set); "
            "solve it instead of lifting"
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wcs = WCS(frame.header, naxis=2)
    return _write_wcs(wcs.to_header(relax=True), cache)


def solve_astap(frame, force=False, astap_exe=ASTAP_EXE, timeout=120, downsample=(2, 0)):
    """Solve and cache a frame's WCS with the local ASTAP solver.

    ASTAP does its own star detection, so it is handed the green plane written to a
    temporary 2-D FITS, with the header pointing and field of view as hints (a large
    robustness and speed win). It solves in place with ``-update``; we then store the
    clean WCS header as our standard sidecar, so this is a drop-in replacement for the
    astrometry.net path.

    Two passes by default: a fast downsampled solve, then a full-resolution retry
    that recovers shallow stacks whose faint stars are lost at ``-z 2``.
    """
    import subprocess
    import tempfile

    cache = wcs_cache_path(frame.path)
    if not force and cache.exists() and _read_wcs(cache).has_celestial:
        return _read_wcs(cache)

    h = frame.header
    ny, nx = frame.shape
    with tempfile.TemporaryDirectory(prefix="astap_") as td:
        tmp = os.path.join(td, "g.fits")
        fits.PrimaryHDU(data=frame.g.astype("float32")).writeto(tmp, overwrite=True)
        cmd = [astap_exe, "-f", tmp, "-update"]
        try:
            cmd += ["-fov", f"{ny * header_pixel_scale(h) / 3600.0:.3f}"]
        except Exception:
            pass  # no optics keywords: let ASTAP solve blind
        try:
            cmd += ["-ra", f"{float(h['RA']) / 15.0:.5f}",
                    "-spd", f"{float(h['DEC']) + 90.0:.5f}"]
        except Exception:
            pass  # no pointing: ditto
        result = None
        for z, radius in zip(downsample, (5, 10)):
            result = subprocess.run(
                cmd + ["-z", str(z), "-r", str(radius)],
                capture_output=True, text=True, timeout=timeout,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                wcs = WCS(fits.getheader(tmp), naxis=2)
            if wcs.has_celestial:
                return _write_wcs(wcs.to_header(relax=True), cache)
        tail = (result.stdout or "").strip().splitlines()[-1:] or [result.stderr[-160:]]
        raise RuntimeError(f"ASTAP no solution (rc={result.returncode}): {tail}")


def solve_nova(frame, x, y, api_key=None, force=False, scale_err=10.0, timeout=300):
    """Solve and cache a frame's WCS from a pixel source list via astrometry.net.

    ``x``/``y`` are 0-indexed source coordinates, brightest first. Sending a source
    list rather than the image is much faster and avoids re-detecting stars we have
    already found.

    Note the service handles only a handful of concurrent jobs per account, so a
    small thread pool (~4) is the practical ceiling; more just produces dropped
    connections.
    """
    cache = wcs_cache_path(frame.path)
    if not force and cache.exists() and _read_wcs(cache).has_celestial:
        return _read_wcs(cache)

    from astroquery.astrometry_net import AstrometryNet

    api_key = api_key or os.environ.get("ASTROMETRY_KEY")
    if not api_key:
        raise ValueError(
            "astrometry.net API key required: set $ASTROMETRY_KEY or pass api_key="
        )

    h = frame.header
    ny, nx = frame.shape
    ast = AstrometryNet()
    ast.api_key = api_key
    # astrometry.net pixel coordinates are 1-indexed; SEP's are 0-indexed.
    wcs_header = ast.solve_from_source_list(
        np.asarray(x) + 1, np.asarray(y) + 1, nx, ny,
        scale_units="arcsecperpix",
        scale_est=header_pixel_scale(h),
        scale_err=scale_err,
        center_ra=float(h["RA"]),
        center_dec=float(h["DEC"]),
        radius=2.0,
        solve_timeout=timeout,
    )
    if not wcs_header:
        raise RuntimeError(f"astrometry.net failed to solve {frame.path}")
    return _write_wcs(wcs_header, cache)


# --- the catalogue-anchored solver ------------------------------------------------------

#: Detections and catalogue sources fed to the bootstrap, brightest first. The catalogue
#: gets more because it covers more sky than the frame does -- see :func:`_bootstrap`.
N_BOOTSTRAP_DETECTIONS = 250
N_BOOTSTRAP_CATALOGUE = 500

#: Fraction of the frame's half-diagonal within which detections are used for the
#: rotation vote. The nominal plate scale is ~1% out, which displaces a star by ``1% * r``
#: -- negligible near the centre and tens of pixels at a corner. Voting on the core keeps
#: the peak sharp; :func:`_best_scale` then recovers the scale from the whole frame.
VOTE_CORE_FRACTION = 0.55

#: Plate-scale factors scanned once the orientation is known. The nominal scale from
#: ``XPIXSZ``/``FOCALLEN`` ran 0.8% high on an S50 and 1.5% high on an S30pro against the
#: solved value, so the grid is centred on 1 and reaches well past both.
#:
#: Ordered outwards from 1.0, and the scan only accepts a *strict* improvement, so a tie
#: keeps the nominal scale. Without that a sparse frame, where many factors pair the same
#: few sources, picks whichever happened to be first -- in practice the end of the grid.
SCALE_FACTORS = 1.0 + np.concatenate([[0.0], np.repeat(np.arange(1, 21), 2)
                                      * np.tile([0.002, -0.002], 20)])

#: Votes at which a parity is accepted without trying the other. Both Seestar models
#: measured positive (see :data:`SEED_PARITY`), so the second is rarely needed and
#: skipping it halves the usual cost.
VOTES_DECISIVE = 40

#: Catalogue sources used in the plate-scale scan, brightest first.
N_SCALE_SCAN = 2000

#: Most detections carried into the solve, brightest first. A deep wide-field co-add
#: yields tens of thousands, and the faint tail is below the reference catalogue's limit
#: -- it cannot pair with anything, and only slows the fit down.
N_DETECT_MAX = 4000

#: Detection threshold for solving, in background sigma. Deliberately far above
#: ``photometry``'s, which is kept low so that the SNR cut defines the sample. See
#: :func:`detect_for_solve`.
SOLVE_THRESH = 5.0

#: Fraction of detections that must pair with a catalogue source for a solve to be
#: believed. See the check at the end of :func:`_refine` for the measurements behind it.
MIN_PAIR_FRACTION = 0.08

#: How far the header pointing may be from the true field centre, in arcmin. This is the
#: translation the bootstrap has to search over.
#:
#: **Not the ~1 arcmin the on-board WCS error suggests.** Measured against ASTAP on real
#: raw subs: 9.0 arcmin on every S50 frame tried, and 9.3 to 40.1 arcmin on the S30pro.
#: The S50's constancy suggests a fixed offset between the reported pointing and the
#: sensor centre rather than mount error, but either way the search has to cover it.
POINTING_SLACK_ARCMIN = 50.0

#: Rotation grid for the bootstrap, in degrees. Alt-Az frames arrive at any orientation,
#: so the whole circle is searched. The step is bounded by the histogram bin: a star at
#: the frame corner moves ``r * step`` under a rotation error, and that has to stay
#: inside a bin or the vote smears across several.
ROTATION_STEP_DEG = 2.0

#: Offset-histogram bin, in pixels. Wide enough to absorb the rotation step above plus
#: the nominal plate scale's error -- measured at 0.8% on the S50 and 1.5% on the S30pro
#: against the solved value, which displaces a corner star by ~16 px.
VOTE_BIN_PX = 32.0

#: Fewest votes in the winning bin for a bootstrap to be believed, whatever the
#: statistics say. A handful of coincidences can look significant against a nearly empty
#: histogram, which is the sparse-field failure mode.
MIN_VOTES_FLOOR = 8

#: How far above the random background the winning bin must stand, in standard
#: deviations of the per-bin Poisson count.
#:
#: A fixed vote count cannot work here: the background scales with the number of pairs
#: and the search area, so a threshold tuned on a 1920x1080 sub (~7 per bin) rejects
#: every solve on a small sparse frame (~0.3 per bin), and one tuned on the small frame
#: accepts noise on the large one.
VOTE_SIGMA = 8.0

#: Sign of ``CDELT1`` to try when seeding, in order.
#:
#: A Seestar's solved CD matrix has a **positive** determinant: the image is mirrored
#: with respect to the usual north-up/east-left convention, so right ascension increases
#: with *x*. Both bundled layouts agree on this, hence ``+1`` first.
#:
#: The other parity is tried as a fallback rather than assumed away, because getting it
#: wrong is unrecoverable rather than merely slow: a similarity transform cannot express
#: a reflection, so the matcher exhausts every triangle it has and fails, no matter how
#: good the data is. Measured on a bundled stack: 0.3 s to match on the right parity,
#: 7 s to fail on the wrong one.
SEED_PARITY = (1.0, -1.0)


def _seed_wcs(frame, parity=SEED_PARITY[0]):
    """A crude TAN WCS from the header pointing and the nominal plate scale.

    Position is good to ~1 arcmin on a real frame and scale to ~1%; the rotation is
    simply assumed to be zero, which it is not -- Alt-Az field rotation reaches tens of
    degrees. Recovering the rotation is the asterism matcher's job; ``parity`` is the one
    thing it cannot recover for itself. Raises if the header carries no pointing at all,
    since there is then nothing to seed from.
    """
    h = frame.header
    if h.get("RA") is None or h.get("DEC") is None:
        raise RuntimeError(
            f"{frame.path} has no RA/DEC in its header, so there is no pointing to "
            "anchor a local solve on; use solver='astap' for a blind solve"
        )
    ny, nx = frame.shape
    scale = header_pixel_scale(h) / 3600.0
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [nx / 2.0, ny / 2.0]
    wcs.wcs.crval = [float(h["RA"]), float(h["DEC"])]
    wcs.wcs.cdelt = [parity * scale, scale]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return wcs


def _header_wcs(frame):
    """The frame's own header WCS, or ``None`` if it has none.

    Deliberately *not* :func:`lift`: this solution is not trusted, only used as a
    starting point that the fit then replaces. Its ~1 arcmin error is irrelevant to a
    matcher working at a few arcsec, and its rotation is worth having.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wcs = WCS(frame.header, naxis=2)
    return wcs if wcs.has_celestial else None


def _pair_up(x, y, cx, cy, tol_px):
    """Nearest-neighbour pairing between detections and predicted positions.

    Returns ``(det_index, cat_index)`` for pairs within ``tol_px``. Where several
    catalogue sources fall on one detection -- routine in a crowded field at 2.4
    arcsec/px -- only the closest is kept, so a blend contributes one pair rather than
    dragging the fit towards the midpoint of a pile-up.
    """
    from scipy.spatial import cKDTree

    tree = cKDTree(np.column_stack([x, y]))
    dist, det = tree.query(np.column_stack([cx, cy]), distance_upper_bound=tol_px)
    hit = np.isfinite(dist)
    det, cat, dist = det[hit], np.flatnonzero(hit), dist[hit]

    order = np.argsort(dist)
    det, cat = det[order], cat[order]
    _, first = np.unique(det, return_index=True)
    return det[first], cat[first]


#: Most pairs handed to a single fit. ``fit_wcs_from_points`` runs a least-squares solve
#: over every point, and a deep wide-field stack pairs up tens of thousands -- measured
#: at 191 s for one S30pro co-add, against 6 s once capped. A TAN fit has six free
#: parameters; several hundred well-spread pairs determine them as well as ten thousand.
N_FIT_MAX = 600


def _fit(x, y, ra, dec, sip_degree=None):
    """Least-squares TAN fit to matched pixel/sky pairs.

    Subsampled with a stride rather than a slice when there are too many, so the pairs
    used stay spread over the frame instead of clustering wherever the bright stars are.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.wcs.utils import fit_wcs_from_points

    x, y = np.asarray(x), np.asarray(y)
    ra, dec = np.asarray(ra), np.asarray(dec)
    if len(x) > N_FIT_MAX:
        keep = slice(None, None, int(np.ceil(len(x) / N_FIT_MAX)))
        x, y, ra, dec = x[keep], y[keep], ra[keep], dec[keep]

    return fit_wcs_from_points(
        (x, y), SkyCoord(ra * u.deg, dec * u.deg),
        projection="TAN", sip_degree=sip_degree,
    )


def _tolerance_ladder(start_arcsec, target_arcsec):
    """Pairing radii stepping down from ``start`` to ``target``, then tightening."""
    ladder, tol = [], float(start_arcsec)
    while tol > target_arcsec:
        ladder.append(tol)
        tol /= 3.0
    return ladder + [target_arcsec, target_arcsec / 2.0, target_arcsec / 3.0]


def _refine(wcs, x, y, ra, dec, scale, tol_arcsec, min_match, sip_degree,
            start_arcsec=None):
    """Iterate match-and-fit from an approximate WCS. ``None`` if it never converges.

    The tolerance tightens each pass: the first is wide enough to survive the seed's
    error, and later ones drop the false pairs that width let in. Returns the last WCS
    that had enough matches, so a final pass that over-tightens does not throw away a
    good solution.

    ``start_arcsec`` sets how wide the first pass is, and has to reflect where the seed
    came from. A bootstrap only locates the field to a histogram bin, and the nominal
    plate scale is ~1% out, which alone displaces a corner star by tens of pixels --
    starting at the target tolerance would pair almost nothing.
    """
    for tol in _tolerance_ladder(start_arcsec or tol_arcsec, tol_arcsec):
        cx, cy = wcs.world_to_pixel_values(ra, dec)
        det, cat = _pair_up(x, y, cx, cy, tol / scale)
        if len(det) < min_match:
            return None
        wcs = _fit(x[det], y[det], ra[cat], dec[cat], sip_degree=sip_degree)

    # Judge the answer at the tolerance that was asked for, not at whatever the last
    # pass happened to use. Returning the best loose fit instead would hand back a
    # solution that is right at the field centre and degrees out at the corners -- which
    # is exactly what a bad plate scale produces, and it reads as success.
    cx, cy = wcs.world_to_pixel_values(ra, dec)
    det, _cat = _pair_up(x, y, cx, cy, tol_arcsec / scale)

    # An absolute floor is not enough on a rich field: a few thousand detections against
    # a few thousand catalogue sources throw up a handful of coincidences within the
    # tolerance whatever the WCS says, so `min_match` alone accepts nonsense. A *wrong*
    # solve pairs almost nothing proportionally -- measured at 2% on an S30pro sub that
    # came out 65 arcsec off, against 18-55% for every correct solve across both models,
    # including a cloud-affected frame where most detections are spurious.
    return wcs if len(det) >= max(min_match, MIN_PAIR_FRACTION * len(x)) else None


def detect_for_solve(frame, thresh=SOLVE_THRESH):
    """Green-plane detections as ``(x, y)``, brightest first. Positions only.

    Solving needs a source *list*, not photometry, and the difference is not small:
    :func:`photometry.extract_sources` sizes apertures per band, measures curves of
    growth and forced-photometers every detection. On a deep S30pro co-add (3840x2160)
    that measured **186 seconds**, none of which affects where the stars are.

    The threshold is also higher than the photometry's. The measurement path keeps it
    low on purpose, so that the SNR cut and not the detection defines the sample; a
    solver wants the opposite. Detecting at 5 sigma instead of 2 took the same S30pro
    frame from 49 s to 8 s, and the sources it drops are the faint ones that the
    reference catalogue does not reach anyway -- they cannot pair with anything, and
    only offer the matcher more ways to go wrong.
    """
    from . import photometry

    _bkg, _sub, objs = photometry._detect(frame.g, thresh)
    order = np.argsort(np.asarray(objs["flux"]))[::-1][:N_DETECT_MAX]
    return np.asarray(objs["x"])[order], np.asarray(objs["y"])[order]


def _rotate(points, centre, degrees):
    """Rotate an ``(n, 2)`` array of pixel coordinates about ``centre``."""
    angle = np.radians(degrees)
    cos, sin = np.cos(angle), np.sin(angle)
    offset = points - centre
    return np.column_stack([
        offset[:, 0] * cos - offset[:, 1] * sin,
        offset[:, 0] * sin + offset[:, 1] * cos,
    ]) + centre


def _vote(detections, predicted, half_px, bin_px=VOTE_BIN_PX):
    """Most popular offset between two point clouds, and how popular it was.

    Every detection is differenced against every predicted position and the offsets are
    histogrammed. If the two clouds share stars under a pure translation, all those pairs
    land in one bin and everything else spreads out flat. Returns ``(votes, dx, dy)``.

    Hand-rolled binning rather than ``np.histogram2d``, which is several times slower and
    is called a few hundred times per frame.
    """
    n_bins = max(int(2 * half_px / bin_px), 1)
    dx = (detections[:, 0][:, None] - predicted[:, 0][None, :]).ravel()
    dy = (detections[:, 1][:, None] - predicted[:, 1][None, :]).ravel()
    ix = np.floor((dx + half_px) * (1.0 / bin_px))
    iy = np.floor((dy + half_px) * (1.0 / bin_px))
    inside = (ix >= 0) & (ix < n_bins) & (iy >= 0) & (iy < n_bins)
    if not inside.any():
        return 0.0, 0.0, 0.0

    flat = (ix[inside] * n_bins + iy[inside]).astype(np.intp)
    counts = np.bincount(flat, minlength=n_bins * n_bins)
    peak = int(counts.argmax())

    # The offset of the winning *bin* is only good to half a bin, and the next step pairs
    # sources at about that radius -- so on a sparse frame the quantisation alone loses
    # most of the matches. Averaging the pairs that actually voted costs nothing and
    # gives the shift to well under a pixel.
    won = flat == peak
    return (float(counts[peak]),
            float(dx[inside][won].mean()),
            float(dy[inside][won].mean()))


def _best_scale(x, y, seed_xy, centre, shift, tol_px):
    """Plate-scale factor that pairs up the most sources, and that pair count.

    The rotation vote deliberately uses only the core of the frame, where the nominal
    scale being ~1% out barely moves a star. That leaves the scale itself unmeasured, and
    it cannot be left at nominal: 1.5% of an S30pro's 2200-pixel half-diagonal is 33
    pixels, which is further than the typical distance to the *wrong* neighbour, so
    pairing at the edges would silently pick the wrong star.

    Scaling happens about the frame centre, which the shift already places correctly, so
    the two are independent and a 1-D scan is enough.
    """
    from scipy.spatial import cKDTree

    tree = cKDTree(np.column_stack([x, y]))
    # Bounded: a wide field can put tens of thousands of catalogue sources in range, and
    # the scan queries the tree once per source per factor. The brightest few thousand
    # settle a single scalar perfectly well.
    offset = (seed_xy - centre)[:N_SCALE_SCAN]
    best = (0, 1.0)
    for factor in SCALE_FACTORS:
        placed = offset * factor + centre + shift
        hits = int(np.isfinite(
            tree.query(placed, distance_upper_bound=tol_px)[0]
        ).sum())
        if hits > best[0]:
            best = (hits, float(factor))
    return best[1], best[0]


def _vote_threshold(n_detections, n_catalogue, half_px, bin_px=VOTE_BIN_PX):
    """Votes the winning bin must draw to be more than a coincidence.

    The background is flat, so each bin holds ``n_pairs / n_bins`` on average and varies
    like its square root. See :data:`VOTE_SIGMA`.
    """
    n_bins = max(int(2 * half_px / bin_px), 1) ** 2
    expected = n_detections * n_catalogue / n_bins
    return max(MIN_VOTES_FLOOR, expected + VOTE_SIGMA * np.sqrt(expected))


def _bootstrap(frame, x, y, ra, dec, scale, min_match=8):
    """Recover an approximate WCS from the pointing alone, by voting on offsets.

    The frame's orientation is unknown -- Alt-Az field rotation reaches tens of degrees
    -- and the pointing is off by up to :data:`POINTING_SLACK_ARCMIN`. What *is* known is
    the plate scale, to about 1%. So the unknown is a rotation and a shift, and both are
    found by brute force: for each angle on a grid, rotate the catalogue's predicted
    pixel positions and histogram every detection-minus-prediction offset. The angle
    whose histogram has the sharpest peak is the orientation, and the peak is the shift.

    An asterism matcher (astroalign, as :mod:`stacking` uses between raw subs) was the
    obvious choice here and does not work, which is worth recording. Triangle invariants
    discard the known scale, so the matcher has to rediscover it, and it can only match
    the brightest N of each list. The catalogue has to cover the frame's diagonal *plus*
    the pointing error, which for an S50 sub is 6.2 square degrees against a 0.92 square
    degree frame -- so only ~15% of the brightest catalogue sources are in the frame at
    all, the two brightness rankings disagree, and matching fails on real frames however
    it is tuned. Voting has no such problem: contamination adds a flat background to the
    histogram rather than competing with the signal.

    Both parities are tried; see :data:`SEED_PARITY`.
    """
    ny, nx = frame.shape
    centre = np.array([nx / 2.0, ny / 2.0])
    half_diagonal = 0.5 * np.hypot(ny, nx)
    # The pointing error plus the frame's own half-diagonal: the furthest a detection can
    # be from where the seed puts its catalogue counterpart.
    half_px = (POINTING_SLACK_ARCMIN * 60.0 / scale) + half_diagonal

    core = np.hypot(x - centre[0], y - centre[1]) <= VOTE_CORE_FRACTION * half_diagonal
    detections = np.column_stack([x[core], y[core]])[:N_BOOTSTRAP_DETECTIONS]
    if len(detections) < min_match:
        raise RuntimeError(
            f"only {len(detections)} detections near the centre of {frame.path}; "
            "too few to anchor a local solve"
        )
    n_catalogue = min(len(ra), N_BOOTSTRAP_CATALOGUE)
    needed = _vote_threshold(len(detections), n_catalogue, half_px)

    best = (0.0, None, None, None)
    for parity in SEED_PARITY:
        seed_xy = np.column_stack(
            _seed_wcs(frame, parity).world_to_pixel_values(ra, dec)
        )
        predicted = seed_xy[:N_BOOTSTRAP_CATALOGUE]
        for degrees in np.arange(0.0, 360.0, ROTATION_STEP_DEG):
            votes, dx, dy = _vote(detections, _rotate(predicted, centre, degrees),
                                  half_px)
            if votes > best[0]:
                best = (votes, parity, degrees, np.array([dx, dy]))
        if best[0] >= max(needed, VOTES_DECISIVE):
            break

    votes, parity, degrees, shift = best
    if votes < needed:
        raise RuntimeError(
            f"could not match {frame.path} against the catalogue: the best offset drew "
            f"{votes:.0f} votes, below the {needed:.0f} needed to stand clear of the "
            f"background ({len(detections)} core detections, {n_catalogue} of "
            f"{len(ra)} catalogue sources used). The pointing may be wrong by more "
            f"than {POINTING_SLACK_ARCMIN:.0f} arcmin, or the catalogue may not cover it."
        )

    # Orientation is settled; now recover the scale over the full frame, then turn the
    # whole thing back into matched pairs and let the least-squares fit produce the WCS.
    # The transform itself is never used, so its conventions never have to be unpicked.
    seed_xy = np.column_stack(
        _seed_wcs(frame, parity).world_to_pixel_values(ra, dec)
    )
    rotated = _rotate(seed_xy, centre, degrees)
    factor, paired = _best_scale(x, y, rotated, centre, shift, VOTE_BIN_PX / 2.0)
    placed = (rotated - centre) * factor + centre + shift

    det, cat = _pair_up(x, y, placed[:, 0], placed[:, 1], VOTE_BIN_PX / 2.0)
    if len(det) < min_match:
        raise RuntimeError(
            f"could not match {frame.path} against the catalogue: the winning offset "
            f"drew {votes:.0f} votes but only {len(det)} sources paired up "
            f"(scale factor {factor:.4f}, {paired} at best)"
        )
    return _fit(x[det], y[det], ra[cat], dec[cat])


def solve_local(frame, catalogue, x=None, y=None, force=False, thresh=SOLVE_THRESH,
                tol_arcsec=3.0, min_match=8, sip_degree=None):
    """Solve and cache a frame's WCS against a local reference catalogue.

    No index files, no subprocess, no network: the catalogue the pipeline already needs
    for its photometry is also everything the astrometry needs. That makes this the only
    solver with no external dependency at all, and it is faster than ASTAP because it
    starts from the pointing instead of searching for it.

    Two routes to an approximate solution, tried in order:

    1. **The header WCS**, when there is one. It is wrong by ~1 arcmin, which is tens of
       pixels -- far too coarse for photometry, but far finer than the few arcsec a
       nearest-neighbour match needs, so most frames converge straight from here.
    2. **Asterism matching**, when there is no header WCS or it fails to converge. This
       is the case for anything :func:`stacking.stack_frame` produced, which drops the
       reference sub's WCS on purpose.

    Either way the approximate solution is only a starting point: the returned WCS comes
    from a least-squares fit to every matched pair, and is rejected outright if too few
    pairs survive.

    Parameters
    ----------
    frame : SeestarFrame
    catalogue : astropy.table.Table
        Reference sources with ``ra``/``dec`` in degrees, covering the field. Ordering
        does not matter; it is sorted by brightness here.
    x, y : array-like, optional
        Detection pixel coordinates. Extracted from the green plane if not given.
    thresh : float
        Detection threshold in background sigma. Note this is *not*
        ``Project.thresh``, which governs the photometry and is deliberately much
        lower -- see :func:`detect_for_solve`.
    tol_arcsec : float
        Initial pairing radius, tightened over the fit iterations. The default is
        generous enough to absorb the seed's error.
    min_match : int
        Fewest pairs an acceptable solution may rest on.
    sip_degree : int, optional
        Fit SIP distortion of this degree. Off by default -- the Seestar's field is
        1 degree across and the residuals do not call for it.
    """
    cache = wcs_cache_path(frame.path)
    if not force and cache.exists() and _read_wcs(cache).has_celestial:
        return _read_wcs(cache)

    if x is None or y is None:
        x, y = detect_for_solve(frame, thresh=thresh)
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)

    ra, dec, n_all = _field_sources(frame, catalogue)
    scale = pixel_scale(frame)

    wcs = None
    header = _header_wcs(frame)
    if header is not None:
        # Wide enough for the arcmin-scale error a header solution can carry, but not so
        # wide that a header wrong by degrees could fake a fit.
        wcs = _refine(header, x, y, ra, dec, scale, tol_arcsec, min_match, sip_degree,
                      start_arcsec=120.0)
    if wcs is None:
        approx = _bootstrap(frame, x, y, ra, dec, scale, min_match=min_match)
        wcs = _refine(approx, x, y, ra, dec, scale, tol_arcsec, min_match, sip_degree,
                      start_arcsec=VOTE_BIN_PX * scale)
    if wcs is None:
        raise RuntimeError(
            f"local solve of {frame.path} did not converge: fewer than {min_match} "
            f"sources paired within {tol_arcsec} arcsec ({len(x)} detections, "
            f"{len(ra)} of {n_all} catalogue sources near the field)"
        )
    return _write_wcs(wcs.to_header(relax=True), cache)


def _field_sources(frame, catalogue):
    """Catalogue sources plausibly on this frame, brightest first.

    Trimmed to a generous circle around the header pointing -- a project's cached
    catalogue covers the whole dithered area, which is far more sky than one frame sees,
    and every extra source is another false-pair candidate for the matcher. Returns
    ``(ra, dec, n_before)``.
    """
    ny, nx = frame.shape
    scale = pixel_scale(frame)
    # Half the diagonal, plus half a degree of slack for the pointing error and for the
    # seed's unknown rotation swinging the corners around.
    radius = 0.5 * np.hypot(ny, nx) * scale / 3600.0 + 0.5

    ra = np.asarray(catalogue["ra"], dtype=float)
    dec = np.asarray(catalogue["dec"], dtype=float)
    mag = _brightness(catalogue)

    h = frame.header
    if h.get("RA") is not None and h.get("DEC") is not None:
        from .gaiadb import separation_deg

        near = separation_deg(float(h["RA"]), float(h["DEC"]), ra, dec) <= radius
        ra, dec, mag = ra[near], dec[near], mag[near]

    order = np.argsort(mag)
    return ra[order], dec[order], len(catalogue)


def _brightness(catalogue):
    """A magnitude to rank catalogue sources by; all-equal if the table has none.

    Masked and non-finite magnitudes sort last rather than being dropped: a source with
    no catalogue magnitude is still a perfectly good astrometric reference, it just
    should not displace a known-bright one from the asterism list.
    """
    for column in ("phot_g_mean_mag", "v_jkc_mag", "g_mag", "mag"):
        if column in catalogue.colnames:
            values = np.ma.getdata(catalogue[column]).astype(float)
            values = np.where(np.ma.getmaskarray(catalogue[column]), np.inf, values)
            return np.where(np.isfinite(values), values, np.inf)
    return np.zeros(len(catalogue), dtype=float)


def solve_from_sources(frame, x, y, api_key=None, force=False, solver="nova",
                       astap_exe=ASTAP_EXE, catalogue=None):
    """Solve using an already-measured source list where the solver can use one.

    ASTAP detects its own stars, so the source list is simply ignored there; this
    exists so ``Extraction.solve_wcs`` can call one function regardless of solver.
    """
    if solver == "astap":
        return solve_astap(frame, force=force, astap_exe=astap_exe)
    if solver == "nova":
        return solve_nova(frame, x, y, api_key=api_key, force=force)
    if solver == "local":
        return solve_local(frame, _require_catalogue(catalogue), x=x, y=y, force=force)
    raise ValueError(
        f"unknown solver {solver!r} (expected 'astap', 'nova' or 'local')"
    )


def _require_catalogue(catalogue):
    """The catalogue, or an error saying how to supply one."""
    if catalogue is None:
        raise ValueError(
            "solver='local' needs a reference catalogue: pass catalogue=, or let "
            "pipeline.solve_all supply it from Project.catalogue()"
        )
    return catalogue


def solve(frame, solver="astap", api_key=None, force=False, astap_exe=ASTAP_EXE,
          thresh=2.0, catalogue=None):
    """Solve (or load cached) a WCS for a frame. The general entry point.

    ``solver="lift"`` copies a trustworthy header solution instead of solving.
    ``"local"`` and ``"nova"`` both extract the frame first -- one to pair against
    ``catalogue``, the other to build the source list it uploads; ``"astap"`` detects
    its own stars and needs neither.
    """
    if not force:
        cached = load_wcs(frame)
        if cached is not None:
            return cached
    if solver == "lift":
        return lift(frame, force=force)
    if solver == "astap":
        return solve_astap(frame, force=force, astap_exe=astap_exe)
    if solver == "local":
        # `thresh` is not forwarded: it is the photometry's detection threshold, and the
        # solver wants a much higher one. See `detect_for_solve`.
        return solve_local(frame, _require_catalogue(catalogue), force=force)
    if solver == "nova":
        from . import photometry

        return photometry.extract_sources(frame, thresh=thresh).solve_wcs(
            api_key=api_key, force=force, solver="nova"
        )
    raise ValueError(
        f"unknown solver {solver!r} (expected 'astap', 'local', 'nova' or 'lift')"
    )


def match_quality(sources, tol_arcsec=2.0):
    """Summarise cross-match separations: the practical test of a WCS solve.

    Returns ``{"median_arcsec", "p90_arcsec", "matched_frac", "n"}`` over the finite
    ``sep_arcsec`` values in a cross-matched table. A good solve sits near ~0.8
    arcsec median with most sources matched; a median approaching the tolerance, or a
    low matched fraction, means the solve is wrong even though it "succeeded".
    """
    sep = np.asarray(sources["sep_arcsec"], dtype=float)
    sep = sep[np.isfinite(sep)]
    if not len(sep):
        return {"median_arcsec": np.nan, "p90_arcsec": np.nan,
                "matched_frac": 0.0, "n": 0}
    return {
        "median_arcsec": float(np.median(sep)),
        "p90_arcsec": float(np.percentile(sep, 90)),
        "matched_frac": float((sep < tol_arcsec).mean()),
        "n": int(len(sep)),
    }
