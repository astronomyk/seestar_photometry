"""Both FITS layouts must normalise to one in-memory frame, and the two header
dialects must resolve to the same physical quantities."""

import numpy as np
import pytest
from astropy.io import fits

from seestar_photometry import frames

from conftest import (
    NX, NY, make_cube, make_header, write_cube, write_mef,
)


def test_cube_layout_normalises(tmp_path, cube):
    path = write_cube(tmp_path / "a.fit", cube)
    frame = frames.load_frame(path)
    assert frame.layout == "cube"
    assert frame.data.shape == (3, NY, NX)
    assert frame.data.dtype == np.float32
    assert frame.shape == (NY, NX)
    assert frame.star_tab is None


def test_mef_layout_normalises(tmp_path, cube):
    path = write_mef(tmp_path / "b.fits", cube)
    frame = frames.load_frame(path)
    assert frame.layout == "mef"
    assert frame.data.shape == (3, NY, NX)
    assert frame.data.dtype == np.float32
    assert frame.shape == (NY, NX)


def test_both_layouts_give_identical_pixels(tmp_path, cube):
    """The whole point of the normalisation: downstream code can't tell them apart."""
    a = frames.load_frame(write_cube(tmp_path / "a.fit", cube))
    b = frames.load_frame(write_mef(tmp_path / "b.fits", cube))
    np.testing.assert_array_equal(a.data, b.data)


def test_channel_last_cube_is_transposed(tmp_path, cube):
    """A ``(ny, nx, 3)`` export must not be read as three rows of image."""
    path = write_cube(tmp_path / "cl.fit", cube, channel_last=True)
    frame = frames.load_frame(path)
    assert frame.data.shape == (3, NY, NX)
    np.testing.assert_allclose(frame.data, cube, rtol=1e-6)


def test_mef_footprint_is_not_a_science_plane(tmp_path, cube):
    """FOOTPRINT is a 2-D image HDU; treating it as a plane would shift the colours."""
    frame = frames.load_frame(write_mef(tmp_path / "f.fits", cube, with_footprint=True))
    # The blue plane must be blue, not a plane of ones.
    assert frame.b.max() > 100.0
    np.testing.assert_allclose(frame.b, cube[2], rtol=1e-6)


def test_mef_planes_resolve_without_extnames(tmp_path, cube):
    """Fall back to the first three 2-D image HDUs when EXTNAMEs are missing."""
    hdus = [fits.PrimaryHDU(header=make_header("mef"))]
    for plane in cube:
        hdus.append(fits.ImageHDU(data=plane.astype(np.float32)))
    path = tmp_path / "noname.fits"
    fits.HDUList(hdus).writeto(path, overwrite=True)
    frame = frames.load_frame(path)
    np.testing.assert_allclose(frame.data, cube, rtol=1e-6)


def test_star_tab_exposed_for_mef(tmp_path, cube):
    frame = frames.load_frame(write_mef(tmp_path / "s.fits", cube, with_startab=True))
    assert frame.star_tab is not None
    assert "x" in frame.star_tab.colnames


def test_shape_comes_from_data_not_header(tmp_path, cube):
    """The MEF primary HDU has no NAXIS1/2 at all, so nothing may read it from there."""
    path = write_mef(tmp_path / "m.fits", cube)
    header = fits.getheader(path)
    assert header.get("NAXIS1") is None
    assert frames.load_frame(path).shape == (NY, NX)


def test_bands_are_distinguishable(cube_frame, cube):
    """r/g/b accessors must map to axis 0 in R, G, B order."""
    np.testing.assert_allclose(cube_frame.r, cube[0], rtol=1e-6)
    np.testing.assert_allclose(cube_frame.g, cube[1], rtol=1e-6)
    np.testing.assert_allclose(cube_frame.b, cube[2], rtol=1e-6)


# --- header dialects ------------------------------------------------------------------

def test_native_dialect_exposure(tmp_path, cube):
    frame = frames.load_frame(
        write_cube(tmp_path / "n.fit", cube, header=make_header("cube", dialect="native"))
    )
    meta = frames.frame_metadata(frame)
    assert meta["n_exp"] == 39
    assert meta["exptime"] == pytest.approx(10.0)         # per sub-exposure
    assert meta["total_exptime"] == pytest.approx(390.0)  # on-sky total
    assert meta["eqmode"] == 0


def test_crowdsky_dialect_exposure(tmp_path, cube):
    """The trap: EXPTIME here is the *total*, so a naive read gives a 410 s sub."""
    frame = frames.load_frame(
        write_mef(tmp_path / "c.fits", cube,
                  header=make_header("mef", dialect="crowdsky"))
    )
    meta = frames.frame_metadata(frame)
    assert meta["n_exp"] == 41
    assert meta["total_exptime"] == pytest.approx(410.0)
    assert meta["exptime"] == pytest.approx(10.0)  # 410 / 41, not 410
    assert meta["obs_end"] is not None


def test_both_dialects_agree_on_physical_quantities(tmp_path, cube):
    """Two files of the same 41x10s stack must report the same per-sub exposure."""
    native = make_header("cube", dialect="native")
    native["STACKCNT"], native["TOTALEXP"], native["EXPTIME"] = 41, 410.0, 10.0
    a = frames.frame_metadata(
        frames.load_frame(write_cube(tmp_path / "x.fit", cube, header=native))
    )
    b = frames.frame_metadata(
        frames.load_frame(write_mef(tmp_path / "y.fits", cube,
                                    header=make_header("mef", dialect="crowdsky")))
    )
    assert a["n_exp"] == b["n_exp"]
    assert a["exptime"] == pytest.approx(b["exptime"])
    assert a["total_exptime"] == pytest.approx(b["total_exptime"])


def test_sparse_header_degrades_to_nan(tmp_path, cube):
    """Another user's minimal header must flow through, not raise."""
    header = make_header("cube")
    for key in ("STACKCNT", "TOTALEXP", "EXPTIME", "EXPOSURE", "EQMODE"):
        header.remove(key, ignore_missing=True)
    frame = frames.load_frame(write_cube(tmp_path / "sparse.fit", cube, header=header))
    meta = frames.frame_metadata(frame)
    assert np.isnan(meta["total_exptime"])
    assert np.isfinite(meta["airmass"])  # still computable from pointing + site + time


def test_unit_id_and_model(cube_frame):
    assert cube_frame.model == "S50"
    assert frames.unit_id(cube_frame.header) == "8a95aa90"


def test_airmass_is_computed(cube_frame):
    """Seestar headers carry no AIRMASS, so it must be derived and sane."""
    airmass = frames.frame_metadata(cube_frame)["airmass"]
    assert 1.0 <= airmass < 5.0


# --- discovery ------------------------------------------------------------------------

def test_localtree_finds_both_extensions(tmp_path, cube):
    (tmp_path / "sub").mkdir()
    write_cube(tmp_path / "one.fit", cube)
    write_mef(tmp_path / "sub" / "two.fits", cube)
    keys = frames.LocalTree(roots=[tmp_path]).keys()
    assert len(keys) == 2
    assert {p.rsplit(".", 1)[1] for p in keys} == {"fit", "fits"}


def test_localtree_dedupes_overlapping_roots(tmp_path, cube):
    """Nested roots must not measure the same frame twice."""
    (tmp_path / "sub").mkdir()
    write_cube(tmp_path / "sub" / "one.fit", cube)
    keys = frames.LocalTree(roots=[tmp_path, tmp_path / "sub"]).keys()
    assert len(keys) == 1


def test_localtree_curate_hook(tmp_path, cube):
    write_cube(tmp_path / "keep.fit", cube)
    write_cube(tmp_path / "drop.fit", cube)
    keys = frames.LocalTree(
        roots=[tmp_path], curate=lambda p: "keep" in str(p)
    ).keys()
    assert len(keys) == 1 and "keep" in keys[0]


def test_localtree_accepts_a_single_root(tmp_path, cube):
    write_cube(tmp_path / "one.fit", cube)
    assert len(frames.LocalTree(roots=tmp_path).keys()) == 1
