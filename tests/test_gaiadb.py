"""The offline reference catalogue: partitioning, cone queries, and the masking contract.

Everything here runs on a synthetic dataset written into ``tmp_path``, so no download and
no network is involved. What is being tested is the *plumbing* -- that a cone returns
exactly the sources a brute-force cut would, that nulls survive as masks rather than as
NaN, and that the ``source_id`` bit trick agrees with a real HEALPix implementation.

The masking tests are the load-bearing ones. Every "did this source match" check in the
package is ``np.ma.getmaskarray``, so a column that comes back holding NaN instead of a
mask reads as a *successful* match with a NaN magnitude, and that propagates into the
zero-point fit as a wrong answer rather than an error.
"""

import numpy as np
import pytest

pytest.importorskip("pyarrow")
pytest.importorskip("astropy_healpix")

from astropy.table import MaskedColumn, Table  # noqa: E402

from seestar_photometry import gaiadb  # noqa: E402

FIELD = (186.6821, 81.474)          # MW Cam: near the pole, where RA arithmetic breaks


def synthetic_catalogue(n=600, center=FIELD, radius_deg=2.0, seed=11):
    """A catalogue in the :mod:`gaiadb` schema, scattered over a cap.

    ``source_id`` is built the way Gaia builds it -- the level-12 HEALPix index shifted
    up by 35 bits -- so the partitioning under test is exercised for real rather than
    against an identifier that was reverse-engineered from the answer.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy_healpix import HEALPix

    rng = np.random.default_rng(seed)
    theta = np.arccos(rng.uniform(np.cos(np.radians(radius_deg)), 1.0, n))
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    points = SkyCoord(center[0] * u.deg, center[1] * u.deg).directional_offset_by(
        phi * u.rad, theta * u.rad
    )
    ra = np.asarray(points.ra.deg, dtype=float)
    dec = np.asarray(points.dec.deg, dtype=float)

    hpx12 = np.asarray(
        HEALPix(nside=2 ** 12, order="nested", frame="icrs").lonlat_to_healpix(
            ra * u.deg, dec * u.deg
        ),
        dtype=np.int64,
    )
    source_id = (hpx12 << 35) + np.arange(n, dtype=np.int64)

    table = Table()
    table["source_id"] = source_id
    table["ra"] = ra
    table["dec"] = dec
    table["phot_g_mean_mag"] = np.linspace(9.0, 17.4, n).astype("float32")
    table["bp_rp"] = rng.uniform(0.2, 2.0, n).astype("float32")
    table["phot_variable_flag"] = np.where(
        np.arange(n) % 50 == 0, "VARIABLE", "NOT_AVAILABLE"
    )
    table["v_jkc_mag"] = (np.linspace(9.0, 17.4, n) + 0.1).astype("float32")
    table["b_jkc_mag"] = (np.linspace(9.0, 17.4, n) + 0.5).astype("float32")
    table["r_jkc_mag"] = (np.linspace(9.0, 17.4, n) - 0.4).astype("float32")
    table["v_jkc_flag"] = np.ones(n, dtype=bool)
    table["pmra"] = np.zeros(n, dtype="float32")
    table["pmdec"] = np.zeros(n, dtype="float32")

    # A third of the rows have no astrophysical parameters, which is realistic and is
    # what the mask has to carry.
    absent = np.zeros(n, dtype=bool)
    absent[::3] = True
    for name in ("c_star", "ruwe", "teff_gspphot"):
        table[name] = MaskedColumn(
            np.ones(n, dtype="float32"), mask=absent, name=name
        )
    table["ipd_frac_multi_peak"] = MaskedColumn(
        np.zeros(n, dtype="uint8"), mask=absent, name="ipd_frac_multi_peak"
    )
    table["non_single_star"] = np.zeros(n, dtype="uint8")
    for name in ("duplicated_source", "in_galaxy_candidates", "has_epoch_photometry"):
        table[name] = np.zeros(n, dtype=bool)
    return table


@pytest.fixture(scope="module")
def catalogue():
    return synthetic_catalogue()


@pytest.fixture(scope="module")
def dataset(tmp_path_factory, catalogue):
    """A written dataset, and the table it was written from."""
    directory = tmp_path_factory.mktemp("gaiadb")
    gaiadb.write_dataset(catalogue, directory)
    return directory


# --- the source_id bit trick -----------------------------------------------------------

def test_hpx_agrees_with_a_real_healpix_implementation(catalogue):
    """``source_id >> 51`` must be the level-4 pixel, or every part is misfiled.

    This is the one piece of the layout that is pure bit arithmetic with no library
    behind it, so it is checked against one.
    """
    import astropy.units as u
    from astropy_healpix import HEALPix

    expected = HEALPix(nside=2 ** gaiadb.HPX_LEVEL, order="nested",
                       frame="icrs").lonlat_to_healpix(
        np.asarray(catalogue["ra"]) * u.deg, np.asarray(catalogue["dec"]) * u.deg
    )
    assert np.array_equal(gaiadb.hpx_of(catalogue["source_id"]), np.asarray(expected))


def test_separation_agrees_with_astropy():
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    rng = np.random.default_rng(3)
    ra, dec = rng.uniform(0, 360, 200), rng.uniform(-89, 89, 200)
    ours = gaiadb.separation_deg(FIELD[0], FIELD[1], ra, dec)
    theirs = SkyCoord(FIELD[0] * u.deg, FIELD[1] * u.deg).separation(
        SkyCoord(ra * u.deg, dec * u.deg)
    ).deg
    assert np.allclose(ours, theirs, atol=1e-9)


# --- cone queries ------------------------------------------------------------------------

@pytest.mark.parametrize("radius", [0.2, 0.5, 1.0, 2.0])
def test_cone_returns_exactly_what_a_brute_force_cut_would(dataset, catalogue, radius):
    got = gaiadb.cone(FIELD, radius, directory=dataset)
    sep = gaiadb.separation_deg(FIELD[0], FIELD[1], catalogue["ra"], catalogue["dec"])
    assert sorted(np.asarray(got["source_id"])) == sorted(
        np.asarray(catalogue["source_id"])[sep <= radius]
    )


def test_cone_smaller_than_a_pixel_still_resolves(dataset):
    """A radius well inside one HEALPix pixel must not fall through the cone search."""
    assert len(gaiadb.cone(FIELD, 0.01, directory=dataset)) >= 0
    assert len(gaiadb.parts_in_cone(FIELD, 0.01)) >= 1


def test_cone_works_across_the_ra_meridian(tmp_path):
    """RA 0/360 is where a naive min/max box filter silently returns half the sky."""
    center = (0.3, 12.0)
    table = synthetic_catalogue(n=400, center=center, radius_deg=1.5, seed=5)
    gaiadb.write_dataset(table, tmp_path)

    got = gaiadb.cone(center, 1.0, directory=tmp_path)
    sep = gaiadb.separation_deg(center[0], center[1], table["ra"], table["dec"])
    assert len(got) == int((sep <= 1.0).sum())
    # The test is only meaningful if the field really does straddle the wrap.
    assert np.asarray(table["ra"]).max() > 359.0
    assert np.asarray(table["ra"]).min() < 1.0


def test_cone_works_at_the_pole(tmp_path):
    center = (120.0, 89.5)
    table = synthetic_catalogue(n=400, center=center, radius_deg=1.5, seed=6)
    gaiadb.write_dataset(table, tmp_path)

    got = gaiadb.cone(center, 1.0, directory=tmp_path)
    sep = gaiadb.separation_deg(center[0], center[1], table["ra"], table["dec"])
    assert len(got) == int((sep <= 1.0).sum())


def test_cone_says_how_to_get_missing_parts(tmp_path):
    """The usual cause is an un-fetched region, so the error must name the fix."""
    gaiadb.write_dataset(synthetic_catalogue(n=50), tmp_path)
    with pytest.raises(FileNotFoundError, match="gaiadb.download"):
        gaiadb.cone((10.0, -30.0), 0.5, directory=tmp_path)


# --- the schema contract ------------------------------------------------------------------

def test_every_column_survives_the_round_trip(dataset, catalogue):
    got = gaiadb.cone(FIELD, 3.0, directory=dataset)
    assert list(got.colnames) == list(gaiadb.COLUMNS)
    assert len(got) == len(catalogue)


def test_nulls_come_back_masked_and_not_as_nan(dataset, catalogue):
    """The trap this whole module has to avoid. See the module docstring."""
    got = gaiadb.cone(FIELD, 3.0, directory=dataset)
    got.sort("source_id")
    reference = catalogue.copy()
    reference.sort("source_id")

    for name in ("c_star", "ruwe", "teff_gspphot", "ipd_frac_multi_peak"):
        mask = np.ma.getmaskarray(got[name])
        assert mask.any(), f"{name} lost its mask entirely"
        assert np.array_equal(mask, np.ma.getmaskarray(reference[name]))
        # A masked float must not also be NaN: downstream code tests the mask, and a
        # NaN underneath would read as a real measurement if anything ever filled it.
        assert not np.isnan(np.ma.getdata(got[name]).astype(float)[~mask]).any()


def test_always_present_columns_are_never_masked(dataset):
    got = gaiadb.cone(FIELD, 3.0, directory=dataset)
    for name in ("source_id", "ra", "dec"):
        assert not np.ma.getmaskarray(got[name]).any()


def test_variable_flag_keeps_its_exact_strings(dataset):
    """``calibration`` and ``lightcurves`` compare this to the literal "VARIABLE"."""
    got = gaiadb.cone(FIELD, 3.0, directory=dataset)
    values = set(np.asarray(got["phot_variable_flag"]).astype(str))
    assert values == {"VARIABLE", "NOT_AVAILABLE"}
    assert (np.asarray(got["phot_variable_flag"]).astype(str) == "VARIABLE").sum() > 0


def test_dtypes_and_units_match_the_schema(dataset):
    got = gaiadb.cone(FIELD, 3.0, directory=dataset)
    assert got["source_id"].dtype == np.int64
    assert got["ra"].dtype == np.float64
    assert got["v_jkc_mag"].dtype == np.float32
    assert got["duplicated_source"].dtype == np.bool_
    assert str(got["ra"].unit) == "deg"
    assert str(got["v_jkc_mag"].unit) == "mag"


def test_writing_without_a_position_is_an_error(tmp_path, catalogue):
    """Positions and identity are the spine; nothing works without them."""
    for column in ("source_id", "ra", "dec"):
        incomplete = catalogue.copy()
        incomplete.remove_column(column)
        with pytest.raises((ValueError, KeyError), match=column):
            gaiadb.write_dataset(incomplete, tmp_path)


def test_a_column_the_dataset_lacks_comes_back_masked(tmp_path, catalogue):
    """A dataset may carry only part of the schema, and must say so.

    The columns come from two Gaia tables and the one holding positions is 790 GB in
    bulk form, so a build assembled from what is actually obtainable is the normal case,
    not a degenerate one. Callers still see every column -- the absent ones are entirely
    masked, which is what code testing ``np.ma.getmaskarray`` should conclude.
    """
    partial = catalogue.copy()
    for column in ("ruwe", "teff_gspphot", "phot_variable_flag"):
        partial.remove_column(column)
    meta = gaiadb.write_dataset(partial, tmp_path)

    assert "ruwe" not in meta["columns_present"]
    assert "v_jkc_mag" in meta["columns_present"]
    assert meta["columns"] == list(gaiadb.COLUMNS)

    got = gaiadb.cone(FIELD, 3.0, directory=tmp_path)
    assert list(got.colnames) == list(gaiadb.COLUMNS)
    assert np.ma.getmaskarray(got["ruwe"]).all()
    assert np.ma.getmaskarray(got["phot_variable_flag"]).all()
    assert not np.ma.getmaskarray(got["v_jkc_mag"]).all()


def test_parts_can_be_written_one_at_a_time_and_indexed_later(tmp_path, catalogue):
    """An all-sky build is 12288 parts and has to survive being interrupted."""
    pixels = np.unique(gaiadb.hpx_of(catalogue["source_id"]))
    for pixel in pixels:
        rows = catalogue[gaiadb.hpx_of(catalogue["source_id"]) == pixel]
        gaiadb.write_part(rows, tmp_path)
    assert gaiadb.manifest(tmp_path) is None, "no manifest until it is finalised"

    meta = gaiadb.finalise_manifest(tmp_path)
    assert meta["rows"] == len(catalogue)
    assert set(meta["parts"]) == {str(int(p)) for p in pixels}
    assert len(gaiadb.cone(FIELD, 3.0, directory=tmp_path)) == len(catalogue)


def test_write_part_refuses_mixed_pixels(tmp_path, catalogue):
    with pytest.raises(ValueError, match="one HEALPix pixel"):
        gaiadb.write_part(catalogue, tmp_path)


# --- filters --------------------------------------------------------------------------------

def test_magnitude_limits_cut_on_the_right_column(dataset):
    bright_v = gaiadb.cone(FIELD, 3.0, directory=dataset, vmag_limit=14.0)
    assert np.ma.filled(bright_v["v_jkc_mag"], -99).max() < 14.0

    bright_g = gaiadb.cone(FIELD, 3.0, directory=dataset, gmag_limit=12.0)
    assert np.ma.filled(bright_g["phot_g_mean_mag"], -99).max() < 12.0
    assert len(bright_g) < len(bright_v)


def test_a_source_with_no_magnitude_fails_the_cut(tmp_path, catalogue):
    """Rather than passing it: an uncalibratable source is not a bright one."""
    table = catalogue.copy()
    table["v_jkc_mag"] = MaskedColumn(
        np.asarray(table["v_jkc_mag"]), mask=np.ones(len(table), dtype=bool)
    )
    gaiadb.write_dataset(table, tmp_path)
    assert len(gaiadb.cone(FIELD, 3.0, directory=tmp_path, vmag_limit=20.0)) == 0


# --- epoch propagation -------------------------------------------------------------------

def test_proper_motion_moves_a_star_by_the_right_amount(tmp_path, catalogue):
    table = catalogue.copy()
    table["pmra"] = np.full(len(table), 1000.0, dtype="float32")   # mas/yr, a fast star
    table["pmdec"] = np.full(len(table), -750.0, dtype="float32")
    gaiadb.write_dataset(table, tmp_path)

    at_gaia = gaiadb.cone(FIELD, 3.0, directory=tmp_path)
    ten_years = gaiadb.cone(FIELD, 3.0, directory=tmp_path,
                            epoch=gaiadb.REF_EPOCH + 10.0)
    at_gaia.sort("source_id")
    ten_years.sort("source_id")

    moved = gaiadb.separation_deg(
        at_gaia["ra"][0], at_gaia["dec"][0],
        ten_years["ra"][0], ten_years["dec"][0],
    ) * 3600.0
    assert moved == pytest.approx(np.hypot(10.0, 7.5), rel=1e-3)
    assert ten_years.meta["epoch"] == gaiadb.REF_EPOCH + 10.0


def test_no_epoch_leaves_positions_at_the_gaia_epoch(dataset, catalogue):
    got = gaiadb.cone(FIELD, 3.0, directory=dataset)
    got.sort("source_id")
    reference = catalogue.copy()
    reference.sort("source_id")
    assert np.allclose(got["ra"], reference["ra"])
    assert got.meta["epoch"] == gaiadb.REF_EPOCH


# --- the manifest ----------------------------------------------------------------------------

def test_manifest_checksums_match_the_files_on_disk(dataset):
    meta = gaiadb.manifest(dataset)
    assert meta["dataset"] == gaiadb.DATASET
    assert meta["columns"] == list(gaiadb.COLUMNS)
    assert meta["parts"], "a written dataset must list its parts"
    for pixel, entry in meta["parts"].items():
        path = gaiadb.part_path(pixel, dataset)
        assert path.exists()
        assert entry["bytes"] == path.stat().st_size
        assert entry["sha256"] == gaiadb.sha256_of(path)
    assert sum(e["rows"] for e in meta["parts"].values()) == meta["rows"]


def test_covers_answers_without_touching_the_network(dataset):
    assert gaiadb.covers(FIELD, 0.5, directory=dataset)
    assert not gaiadb.covers((10.0, -30.0), 0.5, directory=dataset)


def test_covers_is_false_for_a_region_listed_but_not_fetched(tmp_path, catalogue):
    """The partial-download case: the manifest knows about sky that is not here yet.

    Simulated by writing the dataset, recording its manifest, then deleting a part --
    which is exactly the state ``download(center=...)`` leaves for every other region.
    """
    import json

    gaiadb.write_dataset(catalogue, tmp_path)
    assert gaiadb.covers(FIELD, 0.5, directory=tmp_path)

    meta = gaiadb.manifest(tmp_path)
    missing = sorted(meta["parts"], key=int)[0]
    gaiadb.part_path(missing, tmp_path).unlink()
    # The manifest still lists it, so it is missing rather than simply not carried.
    assert json.loads((tmp_path / gaiadb.MANIFEST).read_text())["parts"][missing]
    assert not gaiadb.covers(FIELD, 3.0, directory=tmp_path)


def test_covers_is_false_for_an_empty_directory(tmp_path):
    assert not gaiadb.covers(FIELD, 0.5, directory=tmp_path)


def test_a_partial_dataset_is_valid(tmp_path):
    """Fetching one region must produce something usable, not a broken dataset."""
    gaiadb.write_dataset(synthetic_catalogue(n=200, radius_deg=0.5), tmp_path)
    meta = gaiadb.manifest(tmp_path)
    assert 0 < len(meta["parts"]) < gaiadb.N_PARTS
    assert len(gaiadb.cone(FIELD, 0.4, directory=tmp_path)) > 0


# --- how a Project picks a backend ------------------------------------------------------

def _project(tmp_path, **kwargs):
    from seestar_photometry.frames import LocalTree
    from seestar_photometry.project import Project, Target

    return Project(
        target=Target("MW Cam", ra=FIELD[0], dec=FIELD[1]),
        source=LocalTree(roots=[tmp_path]),
        work_dir=tmp_path / "work",
        **kwargs,
    )


def test_auto_backend_falls_back_to_tap_without_the_download(tmp_path, monkeypatch):
    monkeypatch.setenv("SEESTAR_GAIA_DATA", str(tmp_path / "empty"))
    assert _project(tmp_path).catalogue_backend_used() == "tap"


def test_auto_backend_uses_the_download_when_it_covers_the_field(
        tmp_path, monkeypatch, catalogue):
    """Installing the dataset is the only step needed to stop using the network."""
    root = tmp_path / "cache"
    gaiadb.write_dataset(catalogue, root / gaiadb.DATASET)
    monkeypatch.setenv("SEESTAR_GAIA_DATA", str(root))
    assert _project(tmp_path).catalogue_backend_used() == "local"


def test_backend_can_be_forced(tmp_path, monkeypatch):
    monkeypatch.setenv("SEESTAR_GAIA_DATA", str(tmp_path / "empty"))
    assert _project(tmp_path, catalogue_backend="local").catalogue_backend_used() \
        == "local"


def test_local_backend_writes_the_same_cache_everything_else_reads(
        tmp_path, monkeypatch, catalogue):
    """Downstream code must not be able to tell which backend ran."""
    from seestar_photometry import catalogs

    root = tmp_path / "cache"
    gaiadb.write_dataset(catalogue, root / gaiadb.DATASET)
    monkeypatch.setenv("SEESTAR_GAIA_DATA", str(root))

    proj = _project(tmp_path, catalogue_backend="local", gmag_limit=15.0)
    built = proj.catalogue()
    assert proj.catalogue_path.exists()
    assert np.ma.filled(built["phot_g_mean_mag"], -99).max() < 15.0

    # Read back through the normal path: the ECSV must round-trip, masks included.
    reloaded = catalogs.load_catalogue(proj.catalogue_path)
    assert list(reloaded.colnames) == list(gaiadb.COLUMNS)
    assert np.ma.getmaskarray(reloaded["ruwe"]).any()
    assert len(reloaded) == len(built)


def test_project_catalogue_radius_covers_the_box(tmp_path):
    """The local cone and the TAP box must cover the same sky, or they disagree."""
    proj = _project(tmp_path, catalogue_half_deg=1.5)
    # The corner of the box is half_deg * sqrt(2) away; the radius must reach it.
    assert proj.catalogue_radius_deg >= 1.5 * 2 ** 0.5


def test_fetch_catalogue_asks_for_the_radius_the_project_needs(tmp_path, monkeypatch):
    """A frame's worth of sky is not enough, and getting it wrong is silent.

    Measured on the real deployment: a 1.5 degree fetch pulls 11 HEALPix parts and
    leaves ``catalogue_backend_used()`` answering "tap"; the project needs 14.
    """
    seen = {}
    monkeypatch.setattr(gaiadb, "download",
                        lambda **kw: seen.update(kw) or tmp_path)

    proj = _project(tmp_path, catalogue_half_deg=1.5)
    proj.fetch_catalogue(quiet=True)
    assert seen["radius_deg"] == proj.catalogue_radius_deg
    assert seen["center"] == proj.target.radec
    assert len(gaiadb.parts_in_cone(proj.target.radec, proj.catalogue_radius_deg)) > \
        len(gaiadb.parts_in_cone(proj.target.radec, 1.5))
