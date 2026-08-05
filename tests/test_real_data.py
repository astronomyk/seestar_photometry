"""Tests against the bundled **real** Seestar frames.

The synthetic tests prove the algorithms recover injected truth. These prove the package
copes with what the instrument actually writes -- which is where the surprises live: the
`EXPTIME` dialect collision, the `FOOTPRINT` plane, the Bayer mosaic phase, a WCS that has
to survive a cutout.

Everything runs offline: the frames, their solved `.wcs` sidecars and a trimmed Gaia table
all ship with the package.

Expected values are the ones measured when the dataset was built. They are asserted with
tolerances loose enough not to fail on a library upgrade, but tight enough that a real
regression in the photometry moves them out of range.
"""

import numpy as np
import pytest

from seestar_photometry import (
    astrometry, calibration, debayer, examples, photometry, quality,
)

# The example data is downloaded on demand, so a run on an offline machine skips this
# module rather than failing. `download()` is called once here instead of inside the first
# test, so a slow fetch is not attributed to whichever test happened to run first.
def _have_data():
    if examples.is_downloaded():
        return True
    try:
        examples.download(quiet=True)
        return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(
    not _have_data(),
    reason="example data unavailable (offline?) -- see examples.download()",
)


# --- the bundled dataset itself --------------------------------------------------------

def test_every_advertised_file_is_present():
    """A missing data file must fail loudly here, not mysteriously in a user's doc run."""
    for name in ("stack_c17_15min", "stack_c17_30min", "stack_saturated",
                 "crowdsky_mef", "gaia_mwcam"):
        assert examples.path(name).exists(), name
    assert len(examples.raw_sub_paths()) == 5


def test_unknown_example_lists_what_exists():
    with pytest.raises(FileNotFoundError, match="Available:"):
        examples.path("no_such_frame")


def test_dataset_stays_small():
    """Guard the download size: CI runners fetch this on a cold cache.

    Sums the expected files rather than the whole directory -- the build directory also
    holds the packed archive, and a user cache may hold more than one dataset version.
    """
    directory = examples.data_dir()
    total = sum((directory / name).stat().st_size for name in examples._EXPECTED)
    assert total < 25e6, f"example data grew to {total/1e6:.1f} MB"


def test_download_is_idempotent():
    """A second call must not re-fetch, and must return the same directory."""
    first = examples.download(quiet=True)
    second = examples.download(quiet=True)
    assert first == second == examples.data_dir()
    assert examples.is_downloaded()


def test_checksum_is_pinned():
    """An unpinned checksum would let a replaced asset through silently."""
    assert len(examples.DATA_SHA256) == 64
    assert examples.DATA_VERSION in examples.DATA_URL
    assert examples.DATA_URL.startswith("https://")


def test_env_override_redirects_the_cache(monkeypatch, tmp_path):
    """SEESTAR_PHOTOMETRY_DATA is how an offline machine supplies the data."""
    monkeypatch.setenv("SEESTAR_PHOTOMETRY_DATA", str(tmp_path))
    assert examples.data_dir() == tmp_path
    assert not examples.is_downloaded()      # empty dir -> not satisfied


# --- layouts ---------------------------------------------------------------------------

@pytest.mark.parametrize("loader,layout", [
    (examples.stack, "cube"),
    (examples.stack_deep, "cube"),
    (examples.stack_saturated, "cube"),
    (examples.crowdsky, "mef"),
])
def test_layouts_normalise(loader, layout):
    frame = loader()
    assert frame.layout == layout
    assert frame.data.shape == (3, 1000, 1000)
    assert frame.data.dtype == np.float32
    assert frame.model == "S50"
    # Real sky is never negative-mean or empty.
    assert frame.g.mean() > 100


def test_crowdsky_footprint_is_not_read_as_a_science_plane():
    """The bundled MEF still carries FOOTPRINT precisely so this stays tested."""
    frame = examples.crowdsky()
    # A footprint plane is all ones; a real blue plane is sky-level.
    assert frame.b.mean() > 100
    assert frame.star_tab is not None and len(frame.star_tab) > 0


def test_native_and_crowdsky_header_dialects():
    """The EXPTIME collision, on real headers rather than synthetic ones."""
    from seestar_photometry import frames as frames_mod

    native = frames_mod.frame_metadata(examples.stack())
    crowd = frames_mod.frame_metadata(examples.crowdsky())
    # Native: EXPTIME is per sub (20 s), TOTALEXP the on-sky total.
    assert native["exptime"] == pytest.approx(20.0)
    assert native["total_exptime"] == pytest.approx(760.0)
    assert native["n_exp"] == 38
    # CrowdSky: EXPTIME is the total, so the per-sub value must be derived.
    assert crowd["total_exptime"] > crowd["exptime"]
    assert crowd["obs_end"] is not None      # OB-START/OB-END present
    assert native["obs_end"] is None         # native carries no end time


def test_wcs_survived_the_cutout():
    """CRPIX was shifted when the cutout was made; if that were wrong, nothing matches."""
    frame = examples.stack()
    wcs = astrometry.load_wcs(frame)
    assert wcs is not None and wcs.has_celestial
    ra, dec = wcs.all_pix2world(500, 500, 0)
    # The field centre must still land on MW Cam, within the ~0.35 deg cutout.
    assert abs(float(ra) - examples.FIELD_RA) < 1.0
    assert abs(float(dec) - examples.FIELD_DEC) < 0.4


def test_sidecar_path_handles_the_gz_suffix():
    """`frame.fits.gz` must cache to `frame.wcs`, not `frame.fits.wcs`."""
    path = astrometry.wcs_cache_path(examples.path("stack_c17_15min"))
    assert path.name == "stack_c17_15min.wcs"
    assert path.exists()


# --- raw subs, debayer, stacking --------------------------------------------------------

def test_raw_sub_is_debayered_on_load():
    frame = examples.raw_subs(1)[0]
    assert frame.layout == "bayer"
    assert frame.data.shape == (3, 1000, 1000)
    assert frame.bayer is not None and frame.bayer.ndim == 2


def test_bayer_pattern_orientation_is_right():
    """The GRBG phase, checked the way it was originally established.

    The two green sub-lattices must agree closely, while red and blue do not -- that is
    what proves green sits on (0,0)+(1,1) and there is no row flip. A wrong phase would
    swap colour planes and silently ruin every colour term.
    """
    raw = examples.raw_subs(1)[0].bayer
    med = debayer.channel_medians(raw, "GRBG")
    g00 = float(np.median(raw[0::2, 0::2]))
    g11 = float(np.median(raw[1::2, 1::2]))
    assert abs(g00 - g11) < 10, "the two green lattices disagree -- wrong Bayer phase"
    assert abs(med["R"] - med["G"]) > 10
    assert abs(med["B"] - med["G"]) > 50
    # The documented Seestar raw balance: R/G/B ~ 986/965/1103 ADU.
    assert med["B"] > med["R"] > med["G"]


def test_debayer_tracks_the_native_green_samples():
    """The demosaiced green must follow the sampled green closely.

    Not identically: green is sampled on a quincunx, so the bilinear kernel mixes in the
    neighbouring green samples even at a sampled pixel. Measured 3% median deviation with
    0.97 correlation. What would signal a broken demosaic is loss of correlation, or a
    systematic offset in the median.
    """
    raw = examples.raw_subs(1)[0].bayer
    cube = debayer.debayer(raw, "GRBG")
    native = raw[0::2, 0::2].astype(float)
    interpolated = cube[1][0::2, 0::2]
    assert np.corrcoef(native.ravel(), interpolated.ravel())[0, 1] > 0.9
    assert np.median(np.abs(interpolated - native)) < 0.06 * np.median(native)
    # No systematic brightening or dimming of the plane.
    assert np.median(interpolated) == pytest.approx(np.median(native), rel=0.01)


def test_debayer_rejects_a_2d_frame_without_a_pattern():
    from astropy.io import fits

    header = fits.Header()
    assert not debayer.is_bayer(header, np.zeros((10, 10)))
    header["BAYERPAT"] = "GRBG"
    assert debayer.is_bayer(header, np.zeros((10, 10)))
    assert not debayer.is_bayer(header, np.zeros((3, 10, 10)))   # already 3-plane


def test_unknown_bayer_pattern_raises():
    with pytest.raises(ValueError, match="unknown Bayer pattern"):
        debayer.debayer(np.zeros((10, 10)), "XYZW")


@pytest.mark.parametrize("pattern", sorted(debayer.BAYER_PATTERNS))
def test_every_pattern_produces_a_full_cube(pattern):
    cube = debayer.debayer(np.arange(64, dtype=np.uint16).reshape(8, 8), pattern)
    assert cube.shape == (3, 8, 8)
    assert np.isfinite(cube).all()


@pytest.fixture(scope="module")
def stacked():
    """The five bundled raw subs, registered and co-added once for the whole module."""
    pytest.importorskip("astroalign")
    pytest.importorskip("skimage")
    from seestar_photometry import stacking

    return stacking.stack_frame(examples.raw_sub_paths())


class TestStacking:
    """Local co-add of the five bundled raw subs."""

    def test_all_subs_register(self, stacked):
        _frame, report = stacked
        assert report.n_ok == 5 and report.n_bad == 0
        assert report.resid_median_px < 0.5, "registration residual should be sub-pixel"
        assert report.cover_frac > 0.95

    def test_alt_az_field_rotation_is_detected(self, stacked):
        """These are Alt-Az frames, so the field genuinely rotates between subs.

        Finding ~0 rotation would mean the similarity fit collapsed to a pure shift, which
        is the failure mode that smears stars at the frame corners.
        """
        _frame, report = stacked
        assert 0.05 < report.rot_span_deg < 5.0

    def test_stack_reports_native_exposure_metadata(self, stacked):
        """A local stack must read like an on-board one, or its exposure is misread."""
        from seestar_photometry import frames as frames_mod

        frame, report = stacked
        meta = frames_mod.frame_metadata(frame)
        assert meta["n_exp"] == 5
        assert meta["total_exptime"] == pytest.approx(100.0)   # 5 x 20 s
        assert meta["exptime"] == pytest.approx(20.0)          # per sub, not the total

    def test_stacking_beats_a_single_sub(self, stacked):
        """The point of stacking: the noise must actually come down, roughly as sqrt(N)."""
        frame, report = stacked
        one = photometry.extract_sources(examples.raw_subs(1)[0])
        many = photometry.extract_sources(frame)
        ratio = float(one.rms[1]) / float(many.rms[1])
        assert ratio == pytest.approx(np.sqrt(5), rel=0.35)
        assert len(many.band("G")) > 1.4 * len(one.band("G"))

    def test_stack_drops_the_reference_wcs(self, stacked):
        """The reference sub's WCS does not describe the co-add, so it must not survive.

        Leaving it in place would look solved and reintroduce the ~1 arcmin error the
        package re-solves to avoid.
        """
        frame, _report = stacked
        assert "CRVAL1" not in frame.header
        assert not frame.header.get("PLTSOLVD")

    def test_a_single_sub_stacks_trivially(self):
        pytest.importorskip("astroalign")
        from seestar_photometry import stacking

        frame, report = stacking.stack_frame(examples.raw_sub_paths(1))
        assert report.n_ok == 1 and report.cover_frac == pytest.approx(1.0)

    def test_no_subs_raises(self):
        pytest.importorskip("astroalign")
        from seestar_photometry import stacking

        with pytest.raises(ValueError, match="no sub-exposures"):
            stacking.coadd([])


# --- photometry and calibration on real sky ---------------------------------------------

@pytest.fixture(scope="module")
def measured():
    """``{name: (frame, extraction, calibration)}`` for the three solved stacks."""
    gaia = examples.gaia()
    out = {}
    for name, loader in (("stack", examples.stack),
                         ("deep", examples.stack_deep),
                         ("saturated", examples.stack_saturated)):
        frame = loader()
        wcs = astrometry.load_wcs(frame)
        ext = photometry.extract_sources(frame, enclosed=0.90)
        ext.match_gaia(gaia, wcs=wcs)
        cal = calibration.fit_zeropoint(ext.sources, band="G")
        out[name] = (frame, ext, cal)
    return out


def test_gaia_table_has_what_calibration_needs():
    gaia = examples.gaia()
    assert len(gaia) > 300
    for col in ("source_id", "ra", "dec", "v_jkc_mag", "b_jkc_mag", "r_jkc_mag", "bp_rp"):
        assert col in gaia.colnames
    assert np.isfinite(np.asarray(gaia["v_jkc_mag"], dtype=float)).sum() > 100


def test_detection_finds_a_realistic_number_of_stars(measured):
    _frame, ext, _cal = measured["stack"]
    green = ext.band("G")
    assert 150 < len(green) < 1500
    assert np.all(np.asarray(green["snr"], dtype=float) > 5)


def test_psf_is_chromatic_on_real_data(measured):
    """The finding the per-band aperture exists for, on real optics."""
    _frame, ext, _cal = measured["stack"]
    r, g, b = (float(v) for v in ext.fwhm)
    assert 2.0 < g < 8.0
    assert r > g and b > g
    assert 1.0 < r / g < 1.4
    assert 1.0 < b / g < 1.5


def test_aperture_is_sized_per_band(measured):
    _frame, ext, _cal = measured["stack"]
    assert ext.aperture[0] > ext.aperture[1] < ext.aperture[2]


def test_zeropoint_is_recovered(measured):
    """Measured 23.644 when the dataset was built."""
    _frame, _ext, cal = measured["stack"]
    assert cal.zeropoint == pytest.approx(23.64, abs=0.15)
    assert cal.rms < 0.05
    assert cal.n_stars > 30
    assert cal.colour_label == "B-R (JKC)"


def test_two_stacks_of_the_same_unit_agree_on_the_zeropoint(measured):
    """Same telescope, same night, two integrations: the ZP is a property of the unit.

    They agreed to 0.012 mag when built. A drift here means something has become
    exposure- or depth-dependent that should not be.
    """
    _f1, _e1, shallow = measured["stack"]
    _f2, _e2, deep = measured["deep"]
    assert deep.zeropoint == pytest.approx(shallow.zeropoint, abs=0.08)


def test_longer_exposure_goes_deeper(measured):
    """760 s -> 1460 s must deepen the 5-sigma limit; it measured 17.15 -> 17.69."""
    rows = {}
    for name in ("stack", "deep"):
        _frame, ext, cal = measured[name]
        rows[name] = quality.frame_quality(ext, cal)
    assert rows["deep"]["total_exptime"] > rows["stack"]["total_exptime"]
    gain = rows["deep"]["v_lim_5sigma"] - rows["stack"]["v_lim_5sigma"]
    assert 0.15 < gain < 1.0, f"depth gain {gain:.2f} mag is not plausible"


def test_saturation_limit_is_measured(measured):
    """The saturated example must actually yield a bright limit, not nan.

    Its cutout is deliberately centred on the clipping star -- see
    tools/build_example_data.py.
    """
    frame, ext, cal = measured["saturated"]
    row = quality.frame_quality(ext, cal)
    assert frame.g.max() > 0.9 * 65535
    assert np.isfinite(row["v_sat"])
    assert 6.0 < row["v_sat"] < 11.0
    assert row["v_sat"] < row["v_lim_5sigma"]      # bright limit is brighter than faint


def test_clean_frames_do_not_report_a_saturation_limit(measured):
    """nan is the correct answer when nothing in the field clips."""
    _frame, ext, cal = measured["stack"]
    assert np.isnan(quality.frame_quality(ext, cal)["v_sat"])


def test_astrometry_is_sub_arcsecond(measured):
    """The practical test of the solve; measured 0.51 arcsec median."""
    _frame, ext, _cal = measured["stack"]
    mq = astrometry.match_quality(ext.band("G"))
    assert mq["median_arcsec"] < 1.5
    assert mq["matched_frac"] > 0.3


def test_frame_table_row_is_complete_on_real_data(measured):
    from test_calibration import REQUIRED_COLUMNS

    _frame, ext, cal = measured["stack"]
    row = quality.frame_quality(ext, cal)
    missing = [c for c in REQUIRED_COLUMNS if c not in row]
    assert not missing, f"missing columns: {missing}"
    assert row["bortle"] >= 1
    assert 15.0 < row["sky_sb"] < 23.0


def test_crowdsky_onboard_metrics_can_be_cross_checked():
    """CrowdSky records its own server-side ZP; ours should be in the same ballpark.

    Not equal -- different aperture and reference epoch -- but a wild disagreement would
    mean one of the two is misreading the frame.
    """
    frame = examples.crowdsky()
    onboard = quality.onboard_quality(frame)
    assert onboard, "the bundled CrowdSky frame should carry ZPG/FWHMG"
    ext = photometry.extract_sources(frame)
    ext.match_gaia(examples.gaia(), wcs=astrometry.load_wcs(frame))
    cal = calibration.fit_zeropoint(ext.sources, band="G")
    # Zero points compare directly; measured +0.030 mag apart.
    assert abs(cal.zeropoint - onboard["onboard_zp_G"]) < 0.3

    # FWHM does NOT compare directly: CrowdSky reports arcsec, we report pixels. Convert
    # through the plate scale before comparing, or the two look 2.4x apart.
    scale = astrometry.pixel_scale(frame)
    ours_arcsec = float(ext.fwhm[1]) * scale
    assert ours_arcsec == pytest.approx(onboard["onboard_fwhm_G_arcsec"], rel=0.15)


# --- the project/pipeline path over bundled frames -------------------------------------

def test_pipeline_runs_end_to_end_on_bundled_frames(tmp_path):
    """The batch stages, on real frames, with no network and no solver."""
    from seestar_photometry import lightcurves, pipeline

    proj = examples.project(tmp_path)
    assert len(proj.frames()) == 3            # the three stack_*.fits.gz

    table = pipeline.build_frame_table(proj, workers=1)
    assert len(table) == 3
    assert np.all(np.isfinite(np.asarray(table["zeropoint"], dtype=float)))

    stars, meas = pipeline.build_measurements(proj, workers=1)
    assert len(stars) > 100
    green = meas[np.asarray(meas["band"]) == "G"]
    counts = np.array([int((np.asarray(green["source_id"]) == s).sum())
                       for s in np.unique(green["source_id"])])
    # Forced photometry measures every catalogue source that lands in a frame's footprint,
    # in every such frame -- so a count is 1, 2 or 3, never 0 and never more than 3. It is
    # not always 3 because the three bundled stacks are different pointings (two c17
    # dithers plus a Zcom20 frame), so their footprints only partly overlap.
    assert counts.min() >= 1 and counts.max() == 3
    # The two c17 stacks share most of their sky, so plenty of sources appear in all three.
    assert (counts == 3).sum() > 50

    comps = lightcurves.select_comparisons(
        stars, mag_range=(10.0, 14.0), colour_tol=None, max_sep_arcmin=None
    )
    assert len(comps) > 10
    lc = lightcurves.differential_lightcurve(
        meas, lightcurves.target_id_of(stars), comps, band="G"
    )
    assert len(lc) == 3
    assert np.isfinite(np.asarray(lc["mag"], dtype=float)).all()
