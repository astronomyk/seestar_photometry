"""Aperture sizing and photometry must recover injected truth."""

import numpy as np
import pytest

from seestar_photometry import photometry

from conftest import FWHM, NX, NY


def test_fwhm_recovered_per_band(cube_frame):
    """The measured FWHM must track the injected one, chromatically."""
    for band in ("R", "G", "B"):
        measured = photometry.measure_fwhm(cube_frame, band=band)
        assert measured == pytest.approx(FWHM[band], rel=0.12), band


def test_chromatic_ordering_is_detected(extraction):
    """R and B must come out broader than G, which is what per-band sizing exists for."""
    r, g, b = extraction.fwhm
    assert r > g and b > g
    assert b > r  # injected B/G = 1.20 vs R/G = 1.10


def test_curve_of_growth_is_monotonic_and_saturates(cube_frame):
    cog = photometry.curve_of_growth(cube_frame, band="G")
    frac = np.asarray(cog["flux_frac"], dtype=float)
    # Monotonic while the curve is still rising. Past the plateau the large apertures
    # are dominated by sky noise, so the curve legitimately jitters by a few 0.1% about
    # 1.0 -- requiring monotonicity there would be asserting noise-free sky.
    rising = frac < 0.95
    assert np.all(np.diff(frac[rising]) > 0)
    assert np.all(np.abs(frac[~rising] - 1.0) < 0.03)
    assert frac[-1] == pytest.approx(1.0)     # normalised to the outermost aperture
    assert cog.meta["n_stars"] > 5


def test_aperture_radius_matches_gaussian_expectation(cube_frame):
    """For a Gaussian, the 90% radius is ~1.07 sigma-FWHM; assert the known scaling."""
    cog = photometry.curve_of_growth(cube_frame, band="G")
    radius = photometry.aperture_for_enclosed_flux(cog, 0.90)
    # r_90 = sigma * sqrt(-2 ln 0.1) = 2.146 sigma = 0.911 x FWHM
    assert radius == pytest.approx(0.911 * FWHM["G"], rel=0.15)


def test_aperture_is_chromatic_by_default(extraction):
    """Same enclosed fraction per band means a different radius per band."""
    r, g, b = extraction.aperture
    assert r > g and b > g


def test_explicit_aperture_overrides_sizing(cube_frame):
    ext = photometry.extract_sources(cube_frame, aperture=6.0)
    np.testing.assert_allclose(ext.aperture, [6.0, 6.0, 6.0])


def test_n_fwhm_sizing(cube_frame):
    ext = photometry.extract_sources(cube_frame, n_fwhm=2.0)
    np.testing.assert_allclose(ext.aperture, 2.0 * np.asarray(ext.fwhm), rtol=1e-6)


def test_larger_enclosed_gives_larger_aperture(cube_frame):
    small = photometry.extract_sources(cube_frame, enclosed=0.90)
    large = photometry.extract_sources(cube_frame, enclosed=0.95)
    assert np.all(np.asarray(large.aperture) > np.asarray(small.aperture))


def test_extraction_finds_the_injected_stars(extraction, truth):
    """Every injected star should be detected in green above SNR 5."""
    green = extraction.band("G")
    assert len(green) >= len(truth) * 0.95


def test_detected_positions_match_truth(extraction, truth):
    """Centroids must land on the injected positions to well under a pixel."""
    green = extraction.band("G")
    gx = np.asarray(green["x"], dtype=float)
    gy = np.asarray(green["y"], dtype=float)
    for row in truth:
        d = np.hypot(gx - row["x"], gy - row["y"])
        assert d.min() < 0.5, (row["x"], row["y"], d.min())


def test_sources_carry_all_three_bands(extraction):
    bands = set(np.asarray(extraction.sources["band"]).tolist())
    assert bands == {"R", "G", "B"}


def test_snr_floor_is_enforced(extraction):
    assert np.all(np.asarray(extraction.sources["snr"], dtype=float) > 5.0)


# --- forced photometry ----------------------------------------------------------------

def test_forced_photometry_recovers_relative_fluxes(cube_frame, catalogue, wcs, truth):
    """Fixed-aperture flux must be proportional to injected flux across 4 magnitudes.

    Absolute recovery is not expected -- a 90% aperture keeps 90% of the light by
    construction -- but the *ratio* to truth must be constant, because a
    flux-dependent ratio is exactly the non-linearity that would break a light curve.
    """
    forced = photometry.forced_photometry(
        cube_frame, catalogue["ra"], catalogue["dec"], wcs,
        source_id=catalogue["source_id"],
    )
    green = forced[np.asarray(forced["band"]) == "G"]
    measured = np.asarray(green["flux"], dtype=float)
    injected = np.asarray(truth["flux_G"], dtype=float)
    ratio = measured / injected
    assert np.median(ratio) == pytest.approx(0.90, abs=0.04)
    assert ratio.std() / np.median(ratio) < 0.02


def test_forced_photometry_keeps_off_chip_rows(cube_frame, wcs):
    """A source off the frame must keep its row, flagged, so series never go ragged."""
    ra, dec = wcs.all_pix2world([10.0, -500.0], [10.0, -500.0], 0)
    forced = photometry.forced_photometry(cube_frame, ra, dec, wcs, source_id=[1, 2])
    green = forced[np.asarray(forced["band"]) == "G"]
    assert len(green) == 2
    assert list(np.asarray(green["on_chip"])) == [True, False]


def test_forced_and_detected_use_the_same_aperture(cube_frame, catalogue, wcs):
    """If these ever diverge, forced fluxes leave the zero point's photometric system."""
    ext = photometry.extract_sources(cube_frame, enclosed=0.93)
    forced = photometry.forced_photometry(
        cube_frame, catalogue["ra"], catalogue["dec"], wcs, enclosed=0.93
    )
    np.testing.assert_allclose(ext.aperture, forced.meta["aperture"], rtol=1e-9)


def test_forced_flux_agrees_with_detected_flux(extraction, cube_frame, catalogue, wcs):
    """Measuring at a catalogue position must match detecting the same star."""
    forced = photometry.forced_photometry(
        cube_frame, catalogue["ra"], catalogue["dec"], wcs,
        source_id=catalogue["source_id"],
    )
    fg = forced[np.asarray(forced["band"]) == "G"]
    eg = extraction.band("G")
    fx = np.asarray(fg["x"], dtype=float)
    ex = np.asarray(eg["x"], dtype=float)
    ey = np.asarray(eg["y"], dtype=float)
    fy = np.asarray(fg["y"], dtype=float)
    ratios = []
    for i in range(len(fg)):
        d = np.hypot(ex - fx[i], ey - fy[i])
        j = int(np.argmin(d))
        if d[j] < 1.0:
            ratios.append(float(eg["flux"][j]) / float(fg["flux"][i]))
    assert len(ratios) > 20
    assert np.median(ratios) == pytest.approx(1.0, abs=0.02)


# --- helpers --------------------------------------------------------------------------

def test_instrumental_mag_handles_non_positive_flux():
    """Negative forced fluxes are normal for absent sources; they must become nan."""
    mag = photometry.instrumental_mag([100.0, 0.0, -5.0])
    assert np.isfinite(mag[0])
    assert np.isnan(mag[1]) and np.isnan(mag[2])


def test_mag_error_from_snr():
    assert photometry.mag_error([100.0])[0] == pytest.approx(0.010857)
    assert np.isnan(photometry.mag_error([0.0])[0])


def test_max_in_aperture_flags_bright_peaks(cube_frame, extraction):
    green = extraction.band("G")
    peak = np.asarray(green["max_pix_value"], dtype=float)
    assert np.all(peak > 0)
    assert peak.max() > np.median(peak)


def test_background_fit_recovers_flat_pedestal(extraction):
    from conftest import PEDESTAL

    bg = photometry.fit_background(extraction, band="G")
    assert bg["pedestal"] == pytest.approx(PEDESTAL, rel=0.05)
    assert bg["resid_std"] < 5.0          # a flat sky needs no gradient
    assert bg["coeffs"].shape == (6,)


def test_aperture_correction_is_negative(cube_frame):
    """A partial aperture is fainter than total, so the correction brightens it."""
    cog = photometry.curve_of_growth(cube_frame, band="G")
    radius = photometry.aperture_for_enclosed_flux(cog, 0.90)
    correction = photometry.aperture_correction(cube_frame, radius, cog=cog)
    assert correction == pytest.approx(2.5 * np.log10(0.90), abs=0.02)


def test_extraction_works_on_mef_layout(mef_frame, cube_frame):
    """Photometry must be layout-blind."""
    a = photometry.extract_sources(cube_frame)
    b = photometry.extract_sources(mef_frame)
    np.testing.assert_allclose(a.aperture, b.aperture, rtol=1e-9)
    assert len(a.sources) == len(b.sources)
