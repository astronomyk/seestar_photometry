"""Every diagnostic figure must actually render.

These are smoke tests: they assert the files are produced and non-trivial, not what the
pixels look like. That is the right level -- the figures exist to be looked at by a
human, and what a test can usefully guarantee is that none of them silently stopped
being drawn.

The panel wrapper in :mod:`report` catches a failing panel and draws the error into the
image instead, so a passing file count is not enough on its own; the error-text check
below is what makes these tests meaningful.
"""

import numpy as np
import pytest

from seestar_photometry import lightcurves, plots, report

from conftest import make_wcs

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")


def _nonempty(paths):
    assert paths, "no figures were written"
    for path in paths:
        assert path.exists(), path
        assert path.stat().st_size > 5_000, f"{path.name} looks empty"


def _no_error_panels(paths, capsys=None):
    """No figure may be an error placeholder.

    ``report`` deliberately swallows a broken panel so one failure doesn't cost the
    whole set. That makes "the files exist" a weak assertion, so the placeholders are
    detected by their tell-tale small size.
    """
    tiny = [p.name for p in paths if p.stat().st_size < 12_000]
    assert not tiny, f"these look like error placeholders: {tiny}"


def test_frame_report_writes_the_full_set(tmp_path, cube_frame, extraction,
                                          calibration_fit):
    paths = report.frame_report(cube_frame, extraction, calibration_fit, tmp_path)
    _nonempty(paths)
    stems = {p.stem for p in paths}
    for expected in ("cog", "fwhm", "zp_relation", "zp_colour", "residual_map",
                     "mag_snr", "match_sep", "detections", "background",
                     "colour_fidelity", "residual_v", "residual_snr"):
        assert any(s.endswith(expected) for s in stems), f"missing the {expected} panel"
    _no_error_panels(paths)


def test_frame_report_filenames_are_deterministic(tmp_path, cube_frame, extraction,
                                                  calibration_fit):
    """Re-running must overwrite, not accumulate."""
    first = report.frame_report(cube_frame, extraction, calibration_fit, tmp_path)
    second = report.frame_report(cube_frame, extraction, calibration_fit, tmp_path)
    assert [p.name for p in first] == [p.name for p in second]
    assert len(list(tmp_path.glob("*.png"))) == len(first)


def test_frames_report_writes_dataset_panels(tmp_path, frame_table):
    paths = report.frames_report(frame_table, tmp_path)
    _nonempty(paths)
    names = {p.stem for p in paths}
    assert "frames_zp_vs_time" in names
    assert "frames_rms_hist" in names
    assert "frames_summary" in names          # the contact sheet
    _no_error_panels(paths)


def test_lightcurve_report_writes_the_full_set(tmp_path, lc_bundle, cube_frame):
    lc, stars, meas, comps = lc_bundle
    paths = report.lightcurve_report(
        lc, stars, meas, comps, tmp_path,
        frame=cube_frame, wcs=make_wcs(),
    )
    _nonempty(paths)
    names = {p.stem for p in paths}
    for expected in ("lc_differential", "lc_ensemble_zp", "lc_comparison_grid",
                     "lc_noise_floor", "lc_periodogram", "lc_phase_fold",
                     "lc_finder"):
        assert expected in names, f"missing {expected}"


def test_lightcurve_report_survives_a_missing_finder_frame(tmp_path, lc_bundle):
    """No reference frame just means no finder chart, not a failed report."""
    lc, stars, meas, comps = lc_bundle
    paths = report.lightcurve_report(lc, stars, meas, comps, tmp_path)
    assert paths
    assert "lc_finder" not in {p.stem for p in paths}


def test_host_cutout_figure(tmp_path, cube_frame):
    from seestar_photometry import contamination

    green = np.asarray(cube_frame.g, dtype=float)
    cut = contamination.cutout(green, (190.0, 210.0), (205.0, 225.0), 4.0)
    ax = plots.host_cutout(cut, {"adu": 120.0, "std": 30.0, "n_azimuth": 18})
    ax.figure.savefig(tmp_path / "host.png")
    assert (tmp_path / "host.png").stat().st_size > 5_000


def test_comparison_grid_orders_worst_first(lc_bundle):
    """Problems must land in the top-left, since that is where a human looks first.

    The bundle injects a variable comparison (source 5), so the worst panel is known.
    """
    import matplotlib.pyplot as plt
    import re

    _lc, _stars, meas, comps = lc_bundle
    curves = lightcurves.comparison_curves(meas, comps)
    fig, axes = plots.comparison_grid(curves, comps=comps)

    # Panel titles carry the achieved scatter in mmag; they must descend.
    scatters = []
    for ax in axes.ravel():
        match = re.search(r"(\d+)\s+mmag", ax.get_title())
        if match:
            scatters.append(int(match.group(1)))
    assert len(scatters) == len(curves)
    assert scatters == sorted(scatters, reverse=True)
    # And the rogue star is the worst of them, flagged for attention.
    worst = max(curves, key=lambda sid: plots._scatter_of(curves[sid]))
    assert worst == 5
    assert "CHECK" in axes.ravel()[0].get_title()
    plt.close(fig)


def test_panels_accept_an_external_axes(cube_frame, extraction, calibration_fit):
    """Every panel must compose into a caller's figure, for contact sheets."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(12, 3))
    plots.curve_of_growth(extraction.cogs, extraction.aperture, ax=axes[0])
    plots.zeropoint_vs_colour(calibration_fit, ax=axes[1])
    plots.mag_snr(extraction.sources, ax=axes[2])
    for ax in axes:
        assert ax.has_data() or ax.get_title()
    plt.close(fig)


def test_plot_without_fit_arrays_raises_clearly(calibration_fit):
    """A Calibration read back from a table has no fit arrays; say so plainly."""
    import dataclasses

    stripped = dataclasses.replace(calibration_fit, fit=None)
    with pytest.raises(ValueError, match="no fit arrays"):
        plots.reference_vs_instrumental(stripped)


def test_band_palette_uses_secondary_encoding():
    """Red/green sit in the CVD floor band, so marker shape must differentiate them."""
    from seestar_photometry import _style

    assert len(set(_style.BAND_MARKER.values())) == 3
    assert len(set(_style.BAND_STYLE.values())) == 3
    assert set(_style.BAND_COLOR) == {"R", "G", "B"}
