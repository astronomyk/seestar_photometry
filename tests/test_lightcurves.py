"""Timing, comparison selection and the ensemble differential must behave.

The measurement table here is built synthetically rather than from frames: injecting a
known transparency drift and a known sinusoid is the only way to assert that the
ensemble removes the former and preserves the latter.
"""

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from seestar_photometry import lightcurves

from conftest import FIELD_DEC, FIELD_RA, make_header

N_FRAMES = 40
N_STARS = 12
TARGET_ID = 1
PERIOD = 0.1294           # days -- the MW Cam delta-Scuti period
AMPLITUDE = 0.08          # mag
ZP_DRIFT = 0.35           # mag of transparency variation across the run
NOISE = 0.004             # mag


def synthetic_measurements(seed=3, amplitude=AMPLITUDE, drift=ZP_DRIFT,
                           noise=NOISE, dropout=None, rogue=None):
    """A measurements table with a known signal buried under a known drift.

    ``dropout`` removes a given comparison from the second half of the run, to test
    that a changing ensemble does not shift the zero point. ``rogue`` makes a given
    comparison variable, to test that the per-comparison diagnostic finds it.
    """
    rng = np.random.default_rng(seed)
    times = 2460000.0 + np.linspace(0.0, 0.35, N_FRAMES)
    # Each star's true magnitude; the target sits mid-range.
    true_mag = np.linspace(11.0, 13.0, N_STARS)
    zp = 21.5 - drift * np.linspace(0.0, 1.0, N_FRAMES)

    rows = []
    for f, (t, zp_f) in enumerate(zip(times, zp)):
        for s in range(N_STARS):
            sid = s + 1
            mag = true_mag[s]
            if sid == TARGET_ID:
                mag += amplitude * np.sin(2 * np.pi * t / PERIOD)
            if rogue is not None and sid == rogue:
                mag += 0.10 * np.sin(2 * np.pi * t / 0.07)
            if dropout is not None and sid == dropout and f >= N_FRAMES // 2:
                continue
            m_inst = mag - zp_f + rng.normal(0.0, noise)
            rows.append((sid, f"frame_{f:03d}", "G", t, 1.2,
                         m_inst, noise, 0, True))
    return Table(
        rows=rows,
        names=("source_id", "frame", "band", "bjd_tdb", "airmass",
               "m_inst", "mag_err", "flag", "on_chip"),
    )


def synthetic_stars():
    """A stars table matching :func:`synthetic_measurements`."""
    n = N_STARS
    rng = np.random.default_rng(11)
    return Table({
        "source_id": np.arange(1, n + 1, dtype=np.int64),
        "ra": FIELD_RA + rng.normal(0, 0.05, n),
        "dec": FIELD_DEC + rng.normal(0, 0.05, n),
        "v_jkc_mag": np.linspace(11.0, 13.0, n),
        "bp_rp": np.full(n, 0.8) + rng.normal(0, 0.05, n),
        "phot_variable_flag": np.array(["NOT_AVAILABLE"] * n),
        "sep_target_arcmin": np.concatenate([[0.0], rng.uniform(1, 25, n - 1)]),
        "is_target": np.concatenate([[True], np.zeros(n - 1, dtype=bool)]),
    })


@pytest.fixture
def measurements():
    return synthetic_measurements()


@pytest.fixture
def stars():
    return synthetic_stars()


# --- timing ---------------------------------------------------------------------------

def test_frame_times_uses_the_true_span_when_present():
    """CrowdSky records OB-START/OB-END, so the mid-point must be exact."""
    header = make_header("mef", dialect="crowdsky")
    target = SkyCoord(FIELD_RA * u.deg, FIELD_DEC * u.deg)
    times = lightcurves.frame_times(header, target, total_exptime=410.0)
    assert times["time_source"] == "span"
    # Span is 21:03:55.341 -> 21:14:50.285, so the midpoint is ~21:09:22.8.
    from astropy.time import Time

    expected = Time("2026-07-24T21:09:22.813", format="isot", scale="utc").mjd
    assert times["mjd_mid"] == pytest.approx(expected, abs=1e-6)


def test_frame_times_falls_back_to_half_the_exposure():
    header = make_header("cube", dialect="native")
    target = SkyCoord(FIELD_RA * u.deg, FIELD_DEC * u.deg)
    times = lightcurves.frame_times(header, target, total_exptime=390.0)
    assert times["time_source"] == "exptime"
    assert (times["mjd_mid"] - times["mjd_obs"]) * 86400.0 == pytest.approx(195.0, abs=0.1)


def test_frame_times_produces_a_barycentric_date():
    header = make_header("cube")
    target = SkyCoord(FIELD_RA * u.deg, FIELD_DEC * u.deg)
    times = lightcurves.frame_times(header, target, total_exptime=390.0)
    assert np.isfinite(times["bjd_tdb"])
    # BJD is JD-based, MJD is 2400000.5 lower; the light-travel correction is < 500 s.
    assert abs(times["bjd_tdb"] - (times["mjd_mid"] + 2400000.5)) < 0.01


def test_frame_times_without_a_site_gives_nan_bjd():
    header = make_header("cube")
    header.remove("SITELAT")
    header.remove("SITELONG")
    target = SkyCoord(FIELD_RA * u.deg, FIELD_DEC * u.deg)
    times = lightcurves.frame_times(header, target, total_exptime=390.0)
    assert np.isnan(times["bjd_tdb"])
    assert np.isfinite(times["mjd_mid"])


def test_missing_exptime_does_not_shift_the_epoch():
    header = make_header("cube")
    target = SkyCoord(FIELD_RA * u.deg, FIELD_DEC * u.deg)
    times = lightcurves.frame_times(header, target, total_exptime=np.nan)
    assert times["mjd_mid"] == pytest.approx(times["mjd_obs"])


# --- comparison selection --------------------------------------------------------------

def test_select_comparisons_excludes_the_target(stars):
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    assert TARGET_ID not in np.asarray(comps["source_id"])


def test_select_comparisons_magnitude_cut(stars):
    comps = lightcurves.select_comparisons(stars, dmag=0.3, colour_tol=None)
    target_mag = float(stars["v_jkc_mag"][0])
    assert np.all(np.abs(np.asarray(comps["v_jkc_mag"]) - target_mag) <= 0.3)


def test_select_comparisons_separation_cut(stars):
    comps = lightcurves.select_comparisons(
        stars, dmag=None, colour_tol=None, max_sep_arcmin=10
    )
    assert np.all(np.asarray(comps["sep_target_arcmin"]) <= 10)


def test_select_comparisons_drops_variables(stars):
    stars = stars.copy()
    flag = np.asarray(stars["phot_variable_flag"]).astype("U20")
    flag[3] = "VARIABLE"
    stars["phot_variable_flag"] = flag
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    assert int(stars["source_id"][3]) not in np.asarray(comps["source_id"])


def test_target_id_of_reads_the_flag(stars):
    assert lightcurves.target_id_of(stars) == TARGET_ID


def test_build_stars_flags_an_explicit_target(catalogue):
    stars = lightcurves.build_stars(
        catalogue, catalogue["source_id"], (FIELD_RA, FIELD_DEC),
        target_id=int(catalogue["source_id"][5]),
    )
    flagged = stars[np.asarray(stars["is_target"])]
    assert len(flagged) == 1
    assert int(flagged["source_id"][0]) == int(catalogue["source_id"][5])


def test_build_stars_rejects_an_unmeasured_target(catalogue):
    with pytest.raises(ValueError, match="not measured"):
        lightcurves.build_stars(
            catalogue, catalogue["source_id"], (FIELD_RA, FIELD_DEC),
            target_id=999999,
        )


# --- the ensemble differential ---------------------------------------------------------

def test_differential_removes_the_transparency_drift(measurements, stars):
    """A 0.35 mag drift must not survive into the light curve."""
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    lc = lightcurves.differential_lightcurve(measurements, TARGET_ID, comps)
    assert len(lc) == N_FRAMES
    # The injected sinusoid has amplitude 0.08, so the scatter should be ~0.057
    # (rms of a sine) -- nothing like the 0.35 mag drift that was removed.
    assert np.std(np.asarray(lc["dmag"])) < 0.09
    # And the drift is genuinely gone: no correlation left with frame order.
    trend = np.polyfit(np.arange(len(lc)), np.asarray(lc["dmag"]), 1)[0]
    assert abs(trend * len(lc)) < 0.05


def test_differential_recovers_the_injected_period(measurements, stars):
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    lc = lightcurves.differential_lightcurve(measurements, TARGET_ID, comps)
    pg = lightcurves.periodogram(lc, min_period=0.05, max_period=0.3)
    assert pg["best_period"] == pytest.approx(PERIOD, rel=0.02)
    assert pg["fap"] < 1e-6


def test_differential_recovers_the_injected_amplitude(measurements, stars):
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    lc = lightcurves.differential_lightcurve(measurements, TARGET_ID, comps)
    dmag = np.asarray(lc["dmag"])
    # Peak-to-peak of a sine is 2A; allow for sampling not hitting the extremes.
    assert (dmag.max() - dmag.min()) == pytest.approx(2 * AMPLITUDE, rel=0.20)


def test_constant_star_gives_a_flat_curve(stars):
    """The precision floor test: no injected signal must give scatter at the noise."""
    meas = synthetic_measurements(amplitude=0.0)
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    lc = lightcurves.differential_lightcurve(meas, TARGET_ID, comps)
    assert np.std(np.asarray(lc["dmag"])) < 3 * NOISE


def test_comparison_dropout_does_not_shift_the_zeropoint(stars):
    """Referencing each comparison to its own catalogue magnitude is what buys this.

    A comparison vanishing halfway through must not step the ensemble zero point --
    with a mean-flux-ratio ensemble it would, and the light curve would show a
    discontinuity exactly at the dropout.
    """
    meas = synthetic_measurements(amplitude=0.0, dropout=4)
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    lc = lightcurves.differential_lightcurve(meas, TARGET_ID, comps)
    dmag = np.asarray(lc["dmag"])
    half = len(dmag) // 2
    step = abs(np.mean(dmag[half:]) - np.mean(dmag[:half]))
    assert step < 3 * NOISE, f"ensemble stepped by {step:.4f} mag at the dropout"
    # And the ensemble really did shrink, so the test is exercising the case.
    assert np.asarray(lc["n_comp"])[-1] < np.asarray(lc["n_comp"])[0]


def test_min_comp_drops_underpopulated_frames(measurements, stars):
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)[:2]
    lc = lightcurves.differential_lightcurve(measurements, TARGET_ID, comps, min_comp=3)
    assert len(lc) == 0


def test_lightcurve_reports_the_ensemble_size(measurements, stars):
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    lc = lightcurves.differential_lightcurve(measurements, TARGET_ID, comps)
    assert np.all(np.asarray(lc["n_comp"]) == len(comps))
    assert lc.meta["n_comp_pool"] == len(comps)


def test_flagged_comparisons_shrink_the_ensemble(stars):
    """Flagging two comparisons in one frame must shrink that frame's ensemble only."""
    meas = synthetic_measurements(amplitude=0.0)
    bad = (np.asarray(meas["frame"]) == "frame_010") & np.isin(
        np.asarray(meas["source_id"]), [3, 4]
    )
    meas["flag"][bad] = 1
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    lc = lightcurves.differential_lightcurve(meas, TARGET_ID, comps, max_flag=0)
    assert len(lc) == N_FRAMES                       # the frame is kept
    n_comp = dict(zip(np.asarray(lc["frame"]), np.asarray(lc["n_comp"])))
    assert n_comp["frame_010"] == len(comps) - 2
    assert n_comp["frame_011"] == len(comps)


def test_flagging_the_target_drops_the_frame(stars):
    """If the target itself is unusable there is nothing to plot for that epoch."""
    meas = synthetic_measurements(amplitude=0.0)
    bad = (np.asarray(meas["frame"]) == "frame_010") & (
        np.asarray(meas["source_id"]) == TARGET_ID
    )
    meas["flag"][bad] = 1
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    lc = lightcurves.differential_lightcurve(meas, TARGET_ID, comps, max_flag=0)
    assert len(lc) == N_FRAMES - 1
    assert "frame_010" not in np.asarray(lc["frame"])


# --- the per-comparison diagnostic -----------------------------------------------------

def test_comparison_curves_cover_every_comparison(measurements, stars):
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    curves = lightcurves.comparison_curves(measurements, comps)
    assert set(curves) == set(int(s) for s in comps["source_id"])


def test_comparison_curves_expose_a_rogue_star(stars):
    """The check that earns its keep: a variable comparison must stand out."""
    meas = synthetic_measurements(amplitude=0.0, rogue=5)
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    curves = lightcurves.comparison_curves(meas, comps)
    scatter = {sid: float(np.std(np.asarray(lc["dmag"]))) for sid, lc in curves.items()}
    worst = max(scatter, key=scatter.get)
    assert worst == 5
    others = [v for k, v in scatter.items() if k != 5]
    assert scatter[5] > 3 * np.median(others)


# --- period tools ---------------------------------------------------------------------

def test_phase_fold_wraps_into_unit_interval(measurements, stars):
    comps = lightcurves.select_comparisons(stars, dmag=None, colour_tol=None)
    lc = lightcurves.differential_lightcurve(measurements, TARGET_ID, comps)
    phase, mag = lightcurves.phase_fold(lc, PERIOD)
    assert phase.min() >= 0.0 and phase.max() < 1.0
    assert len(phase) == len(mag) == len(lc)


def test_periodogram_needs_enough_points():
    lc = Table({"time": [1.0, 2.0], "dmag": [0.0, 0.1]})
    with pytest.raises(ValueError, match="at least 5"):
        lightcurves.periodogram(lc)
