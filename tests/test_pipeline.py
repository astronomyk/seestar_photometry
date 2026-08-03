"""The three pipeline stages, end to end on synthetic frames, offline.

Everything is local: frames are written to a temp tree, WCS sidecars are placed
directly rather than solved, and the catalogue cache is pre-populated so no Gaia query
is attempted. That makes this a real end-to-end test of the orchestration -- resume
behaviour, table schemas, error tolerance -- without a network or a plate solver.
"""

import numpy as np
import pytest
from astropy.table import Table

from seestar_photometry import LocalTree, Project, Target, pipeline

from conftest import (
    FIELD_DEC, FIELD_RA, expected_zeropoint, make_catalogue, make_cube,
    write_cube, write_mef, write_wcs_sidecar,
)

N_SYNTH_FRAMES = 4


@pytest.fixture
def dataset(tmp_path, cube, catalogue):
    """A small synthetic project: 3 native frames + 1 CrowdSky frame, all solved."""
    data = tmp_path / "data"
    data.mkdir()
    for i in range(N_SYNTH_FRAMES - 1):
        # A different noise realisation per frame, so the frames are not identical.
        path = write_cube(data / f"synth_{i:02d}.fit", make_cube(seed=100 + i))
        write_wcs_sidecar(path)
    path = write_mef(data / "synth_mef.fits", make_cube(seed=200))
    write_wcs_sidecar(path)

    work = tmp_path / "work"
    proj = Project(
        target=Target("Synthetic Field", ra=FIELD_RA, dec=FIELD_DEC, source_id=1),
        source=LocalTree(roots=[data]),
        work_dir=work,
    )
    proj.ensure_dirs()
    catalogue.write(proj.catalogue_path, format="ascii.ecsv", overwrite=True)
    return proj


# --- project plumbing -----------------------------------------------------------------

def test_project_paths_live_under_work_dir(dataset):
    for path in (dataset.frames_path, dataset.stars_path,
                 dataset.measurements_path, dataset.diagnostics_dir):
        assert dataset.work_dir in path.parents or path.parent == dataset.work_dir


def test_project_finds_every_frame(dataset):
    assert len(dataset.frames()) == N_SYNTH_FRAMES


def test_project_reads_the_cached_catalogue_without_querying(dataset):
    """A pre-populated cache must be read, never re-queried -- the tests are offline."""
    assert len(dataset.catalogue()) > 20


# --- stage 2: the frame table ----------------------------------------------------------

def test_build_frame_table(dataset):
    frames = pipeline.build_frame_table(dataset, workers=1)
    assert len(frames) == N_SYNTH_FRAMES
    assert dataset.frames_path.exists()
    zp = np.asarray(frames["zeropoint"], dtype=float)
    assert np.allclose(zp, expected_zeropoint(0.90), atol=0.03)
    # Both layouts must appear, and be indistinguishable in the science columns.
    assert set(np.asarray(frames["layout"]).tolist()) == {"cube", "mef"}


def test_frame_table_is_idempotent(dataset):
    """A second run must add nothing, and change nothing."""
    first = pipeline.build_frame_table(dataset, workers=1)
    second = pipeline.build_frame_table(dataset, workers=1)
    assert len(second) == len(first)
    np.testing.assert_allclose(
        np.sort(np.asarray(first["zeropoint"], dtype=float)),
        np.sort(np.asarray(second["zeropoint"], dtype=float)),
    )


def test_frame_table_resumes_after_interruption(dataset):
    """Killing a run and restarting must give the same table as running it once.

    Simulated by building with a limit, then building the rest -- which is exactly the
    state a killed run leaves behind, since rows are checkpointed to disk.
    """
    partial = pipeline.build_frame_table(dataset, workers=1, limit=2)
    assert len(partial) == 2
    complete = pipeline.build_frame_table(dataset, workers=1)
    assert len(complete) == N_SYNTH_FRAMES
    assert len(set(np.asarray(complete["path"]).tolist())) == N_SYNTH_FRAMES


def test_frame_table_reports_no_wcs_rather_than_failing(dataset, cube):
    """An unsolved frame must be skipped with a status, not sink the batch."""
    unsolved = write_cube(
        dataset.work_dir.parent / "data" / "unsolved.fit", make_cube(seed=999)
    )
    frames = pipeline.build_frame_table(dataset, workers=1)
    assert len(frames) == N_SYNTH_FRAMES        # the unsolved frame contributed nothing
    assert str(unsolved) not in np.asarray(frames["path"]).tolist()


def test_corrupt_frame_is_recorded_not_fatal(dataset):
    """A truncated download must be logged and skipped."""
    bad = dataset.work_dir.parent / "data" / "truncated.fits"
    bad.write_bytes(b"SIMPLE  =                    T" + b" " * 50)
    frames = pipeline.build_frame_table(dataset, workers=1)
    assert len(frames) == N_SYNTH_FRAMES
    assert dataset.log_path.exists()
    assert "load_error" in dataset.log_path.read_text(encoding="utf-8")


def test_provenance_columns_are_merged(tmp_path, cube, catalogue):
    data = tmp_path / "d"
    data.mkdir()
    write_wcs_sidecar(write_cube(data / "p.fit", cube))
    proj = Project(
        target=Target("F", ra=FIELD_RA, dec=FIELD_DEC),
        source=LocalTree(roots=[data]),
        work_dir=tmp_path / "w",
        provenance=_provenance,
    )
    proj.ensure_dirs()
    catalogue.write(proj.catalogue_path, format="ascii.ecsv", overwrite=True)
    frames = pipeline.build_frame_table(proj, workers=1)
    assert frames["dataset"][0] == "synthetic"


def _provenance(frame):
    """Module-level so it survives pickling to a worker process."""
    return {"dataset": "synthetic"}


# --- stage 3: the measurement tables ---------------------------------------------------

def test_build_measurements(dataset, catalogue):
    stars, meas = pipeline.build_measurements(dataset, workers=1)
    assert dataset.stars_path.exists() and dataset.measurements_path.exists()
    # Every catalogue source, in every band, in every frame.
    assert len(meas) == pytest.approx(len(catalogue) * 3 * N_SYNTH_FRAMES, rel=0.05)
    assert set(np.asarray(meas["band"]).tolist()) == {"R", "G", "B"}
    assert len(stars) == len(catalogue)


def test_measurements_use_the_lightcurve_aperture(dataset):
    """Light curves use the 0.95 aperture, so their zero point differs from the
    frame table's 0.90 one -- if these ever silently coincided, one of the two
    documented aperture choices would not be taking effect."""
    _stars, meas = pipeline.build_measurements(dataset, workers=1)
    zp = np.unique(np.asarray(meas["zeropoint"], dtype=float))
    assert np.allclose(zp, expected_zeropoint(0.95), atol=0.03)


def test_target_is_flagged_in_the_stars_table(dataset):
    stars, _meas = pipeline.build_measurements(dataset, workers=1)
    flagged = stars[np.asarray(stars["is_target"])]
    assert len(flagged) == 1
    assert int(flagged["source_id"][0]) == 1     # the explicit Target.source_id


def test_measurements_carry_times_and_airmass(dataset):
    _stars, meas = pipeline.build_measurements(dataset, workers=1)
    for col in ("mjd_obs", "mjd_mid", "bjd_tdb", "airmass"):
        assert np.isfinite(np.asarray(meas[col], dtype=float)).all(), col


def test_every_source_has_a_row_in_every_frame(dataset):
    """Forced photometry's guarantee: the sample never goes ragged."""
    _stars, meas = pipeline.build_measurements(dataset, workers=1)
    green = meas[np.asarray(meas["band"]) == "G"]
    counts = {}
    for sid in np.unique(np.asarray(green["source_id"])):
        counts[int(sid)] = int((np.asarray(green["source_id"]) == sid).sum())
    assert len(set(counts.values())) == 1, "sources have differing frame counts"


def test_full_chain_produces_a_lightcurve(dataset):
    """The whole point, in one test: frames in, differential light curve out."""
    from seestar_photometry import lightcurves

    pipeline.build_frame_table(dataset, workers=1)
    stars, meas = pipeline.build_measurements(dataset, workers=1)
    comps = lightcurves.select_comparisons(
        stars, dmag=None, colour_tol=None, exclude_variable=False
    )
    lc = lightcurves.differential_lightcurve(
        meas, lightcurves.target_id_of(stars), comps
    )
    assert len(lc) == N_SYNTH_FRAMES
    assert np.isfinite(np.asarray(lc["mag"], dtype=float)).all()
    # A constant synthetic star measured across 4 frames should be flat.
    assert np.std(np.asarray(lc["dmag"], dtype=float)) < 0.02


def test_load_tables_round_trips(dataset):
    pipeline.build_frame_table(dataset, workers=1)
    pipeline.build_measurements(dataset, workers=1)
    frames, stars, meas = pipeline.load_tables(dataset)
    assert all(isinstance(t, Table) for t in (frames, stars, meas))


def test_load_tables_returns_none_before_building(dataset):
    frames, stars, meas = pipeline.load_tables(dataset)
    assert frames is None and stars is None and meas is None
