"""Synthetic Seestar frames, so the whole chain is testable offline in seconds.

Frames are built by injecting Gaussian PSFs of *known* total flux at *known* pixel
positions onto a noisy pedestal, in both FITS layouts and with a deliberately
chromatic PSF (R and B broader than G, as the real optics are). Because the truth is
known, the tests can assert recovery rather than merely "it ran".

The synthetic catalogue is generated *from* the injected fluxes through the exact
calibration relation the pipeline fits::

    V = m_inst + ZP + k * (colour - colour0)

so :func:`calibration.fit_zeropoint` must recover ``ZP_TRUE`` and ``K_TRUE``, and any
change that breaks that relation shows up immediately.
"""

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

# --- the truth the tests assert against ----------------------------------------------

NY, NX = 420, 380          # deliberately non-square, so axis swaps can't hide
FWHM = {"R": 4.4, "G": 4.0, "B": 4.8}   # chromatic: R/G = 1.10, B/G = 1.20
PEDESTAL = 800.0
NOISE = 8.0
ZP_TRUE = 21.5
K_TRUE = 0.12
N_GRID = 7                 # 7x7 star grid
PITCH = 48                 # px between stars: > 2 x the 20 px COG outer radius
MARGIN = 42


def expected_zeropoint(enclosed=0.90):
    """The zero point the pipeline should recover for a given aperture fraction.

    ``ZP_TRUE`` is defined against each star's *total* injected flux, but the standard
    path measures a fixed aperture holding only ``enclosed`` of it. The missing light
    makes every instrumental magnitude fainter by ``-2.5 log10(enclosed)``, and since
    ``ZP = V - m_inst`` the fitted zero point comes out lower by exactly that amount.

    That absorption is the designed behaviour, not an error: because every band is
    sized to the *same* fraction, the offset is identical across bands, so it cancels
    in colours and disappears into the per-frame zero point. This helper is what the
    tests assert against, so the relationship is pinned rather than rediscovered.
    """
    return ZP_TRUE + 2.5 * np.log10(enclosed)

#: A TAN WCS placing the field near the pole, like the real MW Cam data -- which is
#: where naive RA arithmetic breaks, so the tests exercise the awkward case.
FIELD_RA, FIELD_DEC = 186.6821, 81.474
PIXSCALE = 2.39            # arcsec/px, matching an S50


def _gaussian_sigma(fwhm):
    return fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def star_grid():
    """``(x, y)`` pixel positions of the injected stars, well separated and off-edge."""
    coords = [
        (MARGIN + i * PITCH, MARGIN + j * PITCH)
        for i in range(N_GRID) for j in range(N_GRID)
        if MARGIN + i * PITCH < NX - MARGIN and MARGIN + j * PITCH < NY - MARGIN
    ]
    return np.array([c[0] for c in coords], float), np.array([c[1] for c in coords], float)


def truth_table():
    """One row per injected star: position, reference V, colour, and per-band flux.

    Magnitudes span 10-14 (the default fit window) and colours 0.4-1.6, so the
    zero-point fit has real leverage on the colour term rather than fitting noise.
    """
    x, y = star_grid()
    n = len(x)
    rng = np.random.default_rng(20260803)
    v = np.linspace(10.2, 13.8, n)
    colour = 0.4 + 1.2 * rng.random(n)

    # Invert the calibration relation to get the green flux each star must have.
    m_inst_g = v - ZP_TRUE - K_TRUE * (colour - np.median(colour))
    flux_g = 10.0 ** (-m_inst_g / 2.5)
    # R and B carry the colour: B - R = colour, split symmetrically about green so the
    # instrumental colour test has something real to recover.
    flux_r = flux_g * 10.0 ** (0.5 * colour / 2.5)
    flux_b = flux_g * 10.0 ** (-0.5 * colour / 2.5)

    return Table({
        "x": x, "y": y, "v": v, "colour": colour,
        "flux_R": flux_r, "flux_G": flux_g, "flux_B": flux_b,
    })


def make_cube(truth=None, seed=7, chromatic=True):
    """A ``(3, ny, nx)`` float array with the truth table's stars injected."""
    truth = truth_table() if truth is None else truth
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:NY, 0:NX]
    cube = np.empty((3, NY, NX), dtype=np.float32)
    for i, band in enumerate(("R", "G", "B")):
        sigma = _gaussian_sigma(FWHM[band] if chromatic else FWHM["G"])
        plane = PEDESTAL + rng.normal(0.0, NOISE, size=(NY, NX))
        norm = 1.0 / (2.0 * np.pi * sigma ** 2)
        for row in truth:
            # Only integrate near the star; a full-frame exponential per star is slow.
            half = int(np.ceil(6 * sigma))
            x0, x1 = max(int(row["x"]) - half, 0), min(int(row["x"]) + half + 1, NX)
            y0, y1 = max(int(row["y"]) - half, 0), min(int(row["y"]) + half + 1, NY)
            dx = xx[y0:y1, x0:x1] - row["x"]
            dy = yy[y0:y1, x0:x1] - row["y"]
            plane[y0:y1, x0:x1] += row[f"flux_{band}"] * norm * np.exp(
                -(dx ** 2 + dy ** 2) / (2 * sigma ** 2)
            )
        cube[i] = plane
    return cube


def make_header(layout="cube", solved=False, dialect=None):
    """A primary header in one of the two dialects.

    ``dialect`` defaults to ``"native"`` for a cube and ``"crowdsky"`` for a MEF,
    matching how the real files come, but can be forced to test the resolution logic.
    """
    dialect = dialect or ("native" if layout == "cube" else "crowdsky")
    h = fits.Header()
    h["TELESCOP"] = "S50_8a95aa90"
    h["INSTRUME"] = "Seestar S50"
    h["OBJECT"] = "MW Cam"
    h["FILTER"] = "IRCUT"
    h["XPIXSZ"] = 2.9
    h["YPIXSZ"] = 2.9
    h["FOCALLEN"] = 250.0
    h["RA"] = FIELD_RA
    h["DEC"] = FIELD_DEC
    h["SITELAT"] = 48.3157
    h["SITELONG"] = 16.3527
    h["CCD-TEMP"] = 20.3
    h["DATE-OBS"] = "2026-07-24T21:03:55.341"
    if dialect == "native":
        h["EQMODE"] = 0
        h["STACKCNT"] = 39
        h["EXPTIME"] = 10.0
        h["EXPOSURE"] = 10.0
        h["TOTALEXP"] = 390.0
    else:
        # CrowdSky: EXPTIME is the *total*, and the true span is recorded explicitly.
        h["NIMAGES"] = 41
        h["EXPTIME"] = 410.0
        h["EXPOSURE"] = 410.0
        h["OBSTOTAL"] = 654.9
        h["OB-START"] = "2026-07-24 21:03:55.341"
        h["OB-END"] = "2026-07-24 21:14:50.285"
        h["ZPG"] = 22.9
        h["FWHMG"] = 6.6
        h["SKYRMSG"] = 7.1
        h["PIXSCALE"] = 2.3746
    if solved:
        h.update(make_wcs().to_header())
        h["PLTSOLVD"] = True
    return h


def make_wcs():
    """A TAN WCS centred on the synthetic field."""
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.crpix = [NX / 2, NY / 2]
    w.wcs.cdelt = [-PIXSCALE / 3600.0, PIXSCALE / 3600.0]
    w.wcs.crval = [FIELD_RA, FIELD_DEC]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def write_cube(path, cube=None, header=None, channel_last=False):
    """Write the native layout: one 3-D image in the primary HDU."""
    cube = make_cube() if cube is None else cube
    data = np.moveaxis(cube, 0, -1) if channel_last else cube
    fits.PrimaryHDU(
        data=data.astype(np.float32), header=header or make_header("cube")
    ).writeto(path, overwrite=True)
    return path


def write_mef(path, cube=None, header=None, with_footprint=True, with_startab=True):
    """Write the CrowdSky layout: empty primary + named planes (+ FOOTPRINT, STAR-TAB).

    ``FOOTPRINT`` is included by default precisely because it is a trap: it is a 2-D
    image HDU that must never be mistaken for a science plane.
    """
    cube = make_cube() if cube is None else cube
    hdus = [fits.PrimaryHDU(header=header or make_header("mef"))]
    for name, plane in zip(("RED", "GREEN", "BLUE"), cube):
        hdus.append(fits.ImageHDU(data=plane.astype(np.float32), name=name))
    if with_footprint:
        hdus.append(fits.ImageHDU(
            data=np.ones((NY, NX), dtype=np.uint8), name="FOOTPRINT"
        ))
    if with_startab:
        x, y = star_grid()
        hdus.append(fits.BinTableHDU(
            Table({"x": x, "y": y, "flux": np.ones(len(x))}), name="STAR-TAB"
        ))
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


def write_wcs_sidecar(frame_path):
    """Write the ``.wcs`` sidecar the pipeline stages expect, without solving."""
    from seestar_photometry import astrometry

    cache = astrometry.wcs_cache_path(frame_path)
    fits.PrimaryHDU(header=make_wcs().to_header()).writeto(cache, overwrite=True)
    return cache


def make_catalogue(truth=None):
    """A reference catalogue matching the injected stars, in Gaia column names."""
    truth = truth_table() if truth is None else truth
    w = make_wcs()
    ra, dec = w.all_pix2world(
        np.asarray(truth["x"]), np.asarray(truth["y"]), 0
    )
    v = np.asarray(truth["v"], dtype=float)
    colour = np.asarray(truth["colour"], dtype=float)
    return Table({
        "source_id": np.arange(1, len(truth) + 1, dtype=np.int64),
        "ra": ra, "dec": dec,
        "phot_g_mean_mag": v,
        "phot_bp_mean_mag": v + 0.5 * colour,
        "phot_rp_mean_mag": v - 0.5 * colour,
        "bp_rp": colour,
        "phot_variable_flag": np.array(["NOT_AVAILABLE"] * len(truth)),
        "v_jkc_mag": v,
        "b_jkc_mag": v + 0.5 * colour,
        "r_jkc_mag": v - 0.5 * colour,
    })


# --- fixtures -------------------------------------------------------------------------

@pytest.fixture(scope="session")
def truth():
    return truth_table()


@pytest.fixture(scope="session")
def cube(truth):
    return make_cube(truth)


@pytest.fixture(scope="session")
def catalogue(truth):
    return make_catalogue(truth)


@pytest.fixture(scope="session")
def wcs():
    return make_wcs()


@pytest.fixture(scope="session")
def cube_frame(tmp_path_factory, cube):
    """A loaded native-layout frame, with its WCS sidecar already in place."""
    from seestar_photometry import frames

    d = tmp_path_factory.mktemp("cube")
    path = write_cube(d / "synthetic_cube.fit", cube)
    write_wcs_sidecar(path)
    return frames.load_frame(path)


@pytest.fixture(scope="session")
def mef_frame(tmp_path_factory, cube):
    """A loaded CrowdSky-layout frame, with its WCS sidecar already in place."""
    from seestar_photometry import frames

    d = tmp_path_factory.mktemp("mef")
    path = write_mef(d / "synthetic_mef.fits", cube)
    write_wcs_sidecar(path)
    return frames.load_frame(path)


@pytest.fixture(scope="session")
def extraction(cube_frame, catalogue, wcs):
    """A cross-matched extraction of the native frame."""
    from seestar_photometry import photometry

    ext = photometry.extract_sources(cube_frame)
    ext.match_gaia(catalogue, wcs=wcs)
    return ext


@pytest.fixture(scope="session")
def calibration_fit(extraction):
    from seestar_photometry import calibration

    return calibration.fit_zeropoint(extraction.sources, band="G")


@pytest.fixture(scope="session")
def frame_table(extraction, calibration_fit):
    """A synthetic frame table with enough rows and spread to exercise the figures.

    Built by perturbing one real quality row rather than by reducing many frames: the
    dataset-level panels only read the table, so what they need is realistic *columns*
    and a realistic *spread*, not genuinely independent frames.
    """
    from astropy.table import Table

    from seestar_photometry import quality

    base = quality.frame_quality(extraction, calibration_fit)
    rng = np.random.default_rng(42)
    n = 30
    rows = []
    for i in range(n):
        row = dict(base)
        row["frame"] = f"synth_{i:03d}.fit"
        row["path"] = f"/synthetic/synth_{i:03d}.fit"
        row["unit"] = "8a95aa90" if i % 2 else "b1c2d3e4"
        row["date_obs"] = f"2026-07-24T{21 + i // 12:02d}:{(i * 5) % 60:02d}:00"
        row["mjd_mid"] = 60000.0 + i * 0.01
        row["zeropoint"] = base["zeropoint"] + rng.normal(0, 0.08) - 0.004 * i
        row["rms"] = abs(rng.normal(0.04, 0.02))
        row["chi2_red"] = abs(rng.normal(150, 60))
        row["total_exptime"] = float(rng.choice([300.0, 600.0, 900.0, 1800.0]))
        row["fwhm_G"] = base["fwhm_G"] * rng.uniform(0.85, 1.3)
        row["sky_sb"] = rng.uniform(18.5, 21.0)
        row["airmass"] = rng.uniform(1.05, 2.1)
        row["n_cal"] = int(rng.integers(15, 90))
        row["v_lim_5sigma"] = base["v_lim_5sigma"] + rng.normal(0, 0.3)
        rows.append(row)
    return Table(rows)


@pytest.fixture(scope="session")
def lc_bundle():
    """``(lc, stars, measurements, comps)`` for the light-curve figures."""
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from test_lightcurves import synthetic_measurements, synthetic_stars

    from seestar_photometry import lightcurves

    meas = synthetic_measurements(rogue=5)
    stars = synthetic_stars()
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    lc = lightcurves.differential_lightcurve(
        meas, int(stars["source_id"][0]), comps
    )
    return lc, stars, meas, comps
