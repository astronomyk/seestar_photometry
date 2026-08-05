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

#: Catalogue sources and detections fed to the asterism matcher, brightest first.
#: Triangle matching is combinatorial, and the bright end is where the two lists agree --
#: a deep catalogue is mostly stars this frame cannot see. Kept equal to astroalign's
#: ``max_control_points`` below so the truncation happens once, here, and visibly.
N_ASTERISM = 80

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


def _fit(x, y, ra, dec, sip_degree=None):
    """Least-squares TAN fit to matched pixel/sky pairs."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.wcs.utils import fit_wcs_from_points

    return fit_wcs_from_points(
        (np.asarray(x), np.asarray(y)),
        SkyCoord(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg),
        projection="TAN", sip_degree=sip_degree,
    )


def _refine(wcs, x, y, ra, dec, scale, tol_arcsec, min_match, sip_degree):
    """Iterate match-and-fit from an approximate WCS. ``None`` if it never converges.

    The tolerance tightens each pass: the first is wide enough to survive the seed's
    error, and later ones drop the false pairs that width let in. Returns the last WCS
    that had enough matches, so a final pass that over-tightens does not throw away a
    good solution.
    """
    best = None
    for tol in (tol_arcsec, tol_arcsec / 2.0, tol_arcsec / 3.0):
        cx, cy = wcs.world_to_pixel_values(ra, dec)
        det, cat = _pair_up(x, y, cx, cy, tol / scale)
        if len(det) < min_match:
            return best
        wcs = _fit(x[det], y[det], ra[cat], dec[cat], sip_degree=sip_degree)
        best = wcs
    return best


def _bootstrap(frame, x, y, ra, dec):
    """Recover an approximate WCS with no usable orientation, by asterism matching.

    Catalogue sources are projected through :func:`_seed_wcs` into the pixel frame the
    image *would* have at zero rotation, which reduces the problem to matching two point
    clouds differing by a rotation, a small scale error and a shift -- exactly what
    astroalign's triangle matcher solves, and what it already does between raw subs in
    :mod:`stacking`. Both parities are tried; see :data:`SEED_PARITY`.

    The matched pairs go straight into a fit and the fitted transform is discarded, so
    none of astroalign's coordinate conventions have to be reasoned about -- only its
    pairing is used.
    """
    from scipy.spatial import cKDTree

    import astroalign as aa

    detections = np.column_stack([x, y])[:N_ASTERISM]
    failures = []
    for parity in SEED_PARITY:
        seed = _seed_wcs(frame, parity)
        cx, cy = seed.world_to_pixel_values(ra, dec)
        predicted = np.column_stack([cx, cy])
        try:
            _transform, (det_xy, cat_xy) = aa.find_transform(
                detections, predicted[:N_ASTERISM], max_control_points=N_ASTERISM
            )
        except Exception as exc:
            failures.append(f"parity {parity:+.0f}: {type(exc).__name__}: {exc}")
            continue
        _dist, cat = cKDTree(predicted).query(np.asarray(cat_xy))
        return _fit(det_xy[:, 0], det_xy[:, 1], ra[cat], dec[cat])

    # astroalign raises a bare MaxIterError naming neither the frame nor the counts.
    raise RuntimeError(
        f"could not match {frame.path} against the catalogue "
        f"({len(x)} detections, {len(ra)} catalogue sources near the field)\n  "
        + "\n  ".join(failures)
    )


def solve_local(frame, catalogue, x=None, y=None, force=False, thresh=2.0,
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
        from . import photometry

        green = photometry.extract_sources(frame, thresh=thresh).band("G")
        order = np.argsort(np.asarray(green["flux"]))[::-1]  # brightest first
        x, y = np.asarray(green["x"])[order], np.asarray(green["y"])[order]
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)

    ra, dec, n_all = _field_sources(frame, catalogue)
    scale = pixel_scale(frame)

    wcs = None
    header = _header_wcs(frame)
    if header is not None:
        wcs = _refine(header, x, y, ra, dec, scale, tol_arcsec, min_match, sip_degree)
    if wcs is None:
        approx = _bootstrap(frame, x, y, ra, dec)
        wcs = _refine(approx, x, y, ra, dec, scale, tol_arcsec, min_match, sip_degree)
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
        return solve_local(frame, _require_catalogue(catalogue), force=force,
                           thresh=thresh)
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
