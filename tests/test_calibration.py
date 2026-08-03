"""The zero-point fit must recover the injected calibration, and degrade gracefully."""

import numpy as np
import pytest
from astropy.table import MaskedColumn, Table

from seestar_photometry import calibration, quality

from conftest import K_TRUE, ZP_TRUE, expected_zeropoint, make_catalogue, make_wcs


def test_zeropoint_recovered(calibration_fit):
    """The headline test: ZP must come back within a few mmag of the expected value.

    Expected is ``ZP_TRUE`` *minus* the aperture correction -- the 90% aperture's
    missing light is absorbed into the zero point by design. See
    :func:`conftest.expected_zeropoint`.
    """
    assert calibration_fit.zeropoint == pytest.approx(expected_zeropoint(0.90), abs=0.02)


def test_zeropoint_shifts_with_the_aperture_fraction(extraction, cube_frame):
    """A larger aperture keeps more light, so the zero point rises by exactly the
    difference in aperture correction -- and nothing else moves."""
    from seestar_photometry import photometry

    wide = photometry.extract_sources(cube_frame, enclosed=0.95)
    wide.match_gaia(make_catalogue(), wcs=make_wcs())
    cal95 = calibration.fit_zeropoint(wide.sources)
    assert cal95.zeropoint == pytest.approx(expected_zeropoint(0.95), abs=0.02)
    assert cal95.colour_term == pytest.approx(K_TRUE, abs=0.03)


def test_colour_term_recovered(calibration_fit):
    assert calibration_fit.colour_term == pytest.approx(K_TRUE, abs=0.03)


def test_scatter_is_small_for_noise_free_truth(calibration_fit):
    """With a clean synthetic field the fit residuals should be well inside the
    photometric-grade threshold."""
    assert calibration_fit.rms < 0.02


def test_colour0_decorrelates_the_zeropoint(calibration_fit, extraction):
    """ZP is quoted at the median colour, so colour0 must be that median."""
    green = extraction.band("G")
    colour = (np.asarray(green["b_jkc_mag"], dtype=float)
              - np.asarray(green["r_jkc_mag"], dtype=float))
    assert calibration_fit.colour0 == pytest.approx(np.median(colour), abs=0.05)


def test_fit_respects_the_magnitude_window(extraction):
    """Narrowing the window must reduce the star count without moving the zero point."""
    wide = calibration.fit_zeropoint(extraction.sources, mag_range=(10.0, 14.0))
    narrow = calibration.fit_zeropoint(extraction.sources, mag_range=(11.0, 12.5))
    assert narrow.n_stars < wide.n_stars
    assert narrow.zeropoint == pytest.approx(wide.zeropoint, abs=0.05)


def test_colour_cascade_falls_back_to_bp_rp(extraction):
    """Without synthetic JKC colours the fit must fall back to Gaia native BP-RP."""
    sources = extraction.sources.copy()
    sources.remove_columns(["b_jkc_mag", "r_jkc_mag"])
    cal = calibration.fit_zeropoint(sources)
    assert cal.colour_label == "BP-RP"
    assert cal.zeropoint == pytest.approx(expected_zeropoint(), abs=0.05)


def test_fit_degrades_to_a_clipped_mean_without_colour(extraction):
    """No colour at all must still give a usable zero point, with k pinned at zero."""
    sources = extraction.sources.copy()
    sources.remove_columns(["b_jkc_mag", "r_jkc_mag", "bp_rp"])
    cal = calibration.fit_zeropoint(sources)
    assert cal.colour_label == "none"
    assert cal.colour_term == 0.0
    assert np.isnan(cal.colour0)
    # Without the colour term the scatter grows, but the mean must stay right.
    assert cal.zeropoint == pytest.approx(expected_zeropoint(), abs=0.10)


def test_v_column_falls_back_to_v_mag(extraction):
    """A non-Gaia reference catalogue using plain 'v_mag' must work."""
    sources = extraction.sources.copy()
    sources.rename_column("v_jkc_mag", "v_mag")
    cal = calibration.fit_zeropoint(sources)
    assert cal.zeropoint == pytest.approx(expected_zeropoint(), abs=0.02)


def test_missing_v_column_raises(extraction):
    sources = extraction.sources.copy()
    sources.remove_column("v_jkc_mag")
    with pytest.raises(KeyError):
        calibration.fit_zeropoint(sources)


def test_variables_are_excluded(extraction):
    sources = extraction.sources.copy()
    flag = np.asarray(sources["phot_variable_flag"]).astype("U20")
    flag[: len(flag) // 2] = "VARIABLE"
    sources["phot_variable_flag"] = flag
    with_var = calibration.fit_zeropoint(sources, exclude_variable=False)
    without = calibration.fit_zeropoint(sources, exclude_variable=True)
    assert without.n_stars < with_var.n_stars


def test_sigma_clipping_rejects_an_outlier(extraction):
    """One badly wrong star must not move the zero point."""
    sources = extraction.sources.copy()
    flux = np.asarray(sources["flux"], dtype=float)
    green = np.asarray(sources["band"]) == "G"
    idx = np.where(green)[0][:3]
    flux[idx] *= 10.0            # three stars 2.5 mag too bright
    sources["flux"] = flux
    cal = calibration.fit_zeropoint(sources)
    assert cal.zeropoint == pytest.approx(expected_zeropoint(), abs=0.05)


def test_apply_calibration_round_trips(calibration_fit, extraction):
    """Putting instrumental magnitudes back on the V scale must recover reference V."""
    green = extraction.band("G")
    matched = ~np.ma.getmaskarray(green["v_jkc_mag"])
    m_inst = -2.5 * np.log10(np.asarray(green["flux"], dtype=float))
    colour = (np.asarray(green["b_jkc_mag"], dtype=float)
              - np.asarray(green["r_jkc_mag"], dtype=float))
    v = calibration.apply_calibration(m_inst, calibration_fit, colour)
    truth_v = np.asarray(green["v_jkc_mag"], dtype=float)
    residual = (v - truth_v)[matched]
    assert np.nanstd(residual) < 0.02


def test_limiting_mag_gets_fainter_with_lower_noise():
    bright = calibration.limiting_mag(21.5, 100.0)
    faint = calibration.limiting_mag(21.5, 10.0)
    assert faint > bright
    assert faint - bright == pytest.approx(2.5, abs=0.01)  # one dex of noise


def test_limiting_mag_nsigma_scaling():
    five = calibration.limiting_mag(21.5, 10.0, nsigma=5.0)
    hundred = calibration.limiting_mag(21.5, 10.0, nsigma=100.0)
    assert five > hundred
    assert five - hundred == pytest.approx(2.5 * np.log10(20.0), abs=1e-6)


def test_saturation_mag_finds_the_faintest_saturated_star():
    sources = Table({
        "max_pix_value": [65535.0, 65000.0, 1000.0, 500.0],
        "v_jkc_mag": MaskedColumn([8.0, 9.5, 12.0, 13.0], mask=[0, 0, 0, 0]),
    })
    assert calibration.saturation_mag(sources) == pytest.approx(9.5)


def test_saturation_mag_is_nan_when_nothing_saturates():
    sources = Table({
        "max_pix_value": [1000.0, 500.0],
        "v_jkc_mag": MaskedColumn([12.0, 13.0], mask=[0, 0]),
    })
    assert np.isnan(calibration.saturation_mag(sources))


def test_sky_surface_brightness():
    # ZP 21.5, 800 ADU/px over (2.39")^2 -> 21.5 - 2.5 log10(800 / 5.71) = 16.13
    mu = calibration.sky_surface_brightness(800.0, 21.5, 2.39)
    assert mu == pytest.approx(16.13, abs=0.02)
    assert np.isnan(calibration.sky_surface_brightness(-5.0, 21.5, 2.39))


def test_bortle_classes_are_ordered():
    assert calibration.effective_bortle(22.0) == 1
    assert calibration.effective_bortle(20.0) == 5
    assert calibration.effective_bortle(16.0) == 9
    assert np.isnan(calibration.effective_bortle(np.nan))


# --- the frame-table row --------------------------------------------------------------

#: Columns downstream code and the docs rely on. Renaming one silently breaks the
#: figures and the depth model, so the schema is pinned.
REQUIRED_COLUMNS = (
    "path", "frame", "layout", "model", "unit", "telescope", "date_obs",
    "n_exp", "exptime", "total_exptime", "airmass", "pixscale",
    "n_sources", "n_green", "n_matched", "n_cal",
    "zeropoint", "zeropoint_err", "colour_term", "colour_term_err", "colour0",
    "rms", "chi2_red", "v_lim_5sigma", "v_lim_100sigma", "v_sat", "sigma_aper",
    "fwhm_R", "fwhm_G", "fwhm_B", "aperture_R", "aperture_G", "aperture_B",
    "sky_R", "sky_G", "sky_B", "sky_pedestal", "bg_poly", "bg_resid_std",
    "sky_sb", "bortle",
)


def test_frame_quality_schema(extraction, calibration_fit):
    row = quality.frame_quality(extraction, calibration_fit)
    missing = [c for c in REQUIRED_COLUMNS if c not in row]
    assert not missing, f"frame table lost columns: {missing}"


def test_frame_quality_values_are_sane(extraction, calibration_fit):
    row = quality.frame_quality(extraction, calibration_fit)
    assert row["zeropoint"] == pytest.approx(expected_zeropoint(), abs=0.02)
    assert row["v_lim_5sigma"] > row["v_lim_100sigma"]  # 5 sigma reaches fainter
    assert row["fwhm_R"] > row["fwhm_G"] < row["fwhm_B"]
    assert row["n_cal"] <= row["n_matched"] <= row["n_green"]
    assert row["model"] == "S50"
    assert row["unit"] == "8a95aa90"


def test_provenance_hook_adds_columns(extraction, calibration_fit):
    row = quality.frame_quality(
        extraction, calibration_fit,
        provenance=lambda frame: {"dataset": "synthetic", "bin_min": 15.0},
    )
    assert row["dataset"] == "synthetic"
    assert row["bin_min"] == 15.0


def test_onboard_quality_lifted_from_crowdsky_header(mef_frame, cube_frame):
    """A CrowdSky frame carries the server's own metrics; a native one does not."""
    lifted = quality.onboard_quality(mef_frame)
    assert lifted["onboard_zp_G"] == pytest.approx(22.9)
    assert lifted["onboard_fwhm_G"] == pytest.approx(6.6)
    assert quality.onboard_quality(cube_frame) == {}
