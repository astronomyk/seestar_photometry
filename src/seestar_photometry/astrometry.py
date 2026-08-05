"""Per-frame WCS, cached as a ``.wcs`` sidecar next to the FITS file.

The on-board Seestar WCS is **not** photometric-grade: positions are off by tens of
pixels (~1 arcmin), which drops the catalogue match rate to a few percent and makes
forced photometry meaningless. Every frame is therefore re-solved, and the result is
cached as a header-only FITS sidecar so the (expensive) solve happens exactly once
per frame and is shared between projects.

Two solvers, same sidecar, chosen with ``solver=``:

``"astap"`` (default)
    A local, fully offline blind plate solver. Does its own star detection, handles
    the field rotation of Alt-Az frames, and takes ~a second per frame. This is the
    right default: it has no rate limit, no network failure mode, and no API key.

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


def solve_from_sources(frame, x, y, api_key=None, force=False, solver="nova",
                       astap_exe=ASTAP_EXE):
    """Solve using an already-measured source list where the solver can use one.

    ASTAP detects its own stars, so the source list is simply ignored there; this
    exists so ``Extraction.solve_wcs`` can call one function regardless of solver.
    """
    if solver == "astap":
        return solve_astap(frame, force=force, astap_exe=astap_exe)
    if solver == "nova":
        return solve_nova(frame, x, y, api_key=api_key, force=force)
    raise ValueError(f"unknown solver {solver!r} (expected 'astap' or 'nova')")


def solve(frame, solver="astap", api_key=None, force=False, astap_exe=ASTAP_EXE,
          thresh=2.0):
    """Solve (or load cached) a WCS for a frame. The general entry point.

    ``solver="lift"`` copies a trustworthy header solution instead of solving.
    Otherwise, for ``"nova"`` the frame is extracted first to build the source list;
    for ``"astap"`` no extraction is needed.
    """
    if not force:
        cached = load_wcs(frame)
        if cached is not None:
            return cached
    if solver == "lift":
        return lift(frame, force=force)
    if solver == "astap":
        return solve_astap(frame, force=force, astap_exe=astap_exe)
    if solver == "nova":
        from . import photometry

        return photometry.extract_sources(frame, thresh=thresh).solve_wcs(
            api_key=api_key, force=force, solver="nova"
        )
    raise ValueError(f"unknown solver {solver!r} (expected 'astap', 'nova' or 'lift')")


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
