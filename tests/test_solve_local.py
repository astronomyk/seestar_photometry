"""The catalogue-anchored plate solver, on synthetic frames.

Injected stars at known pixel positions and a catalogue derived from the same truth
table, so the solve can be checked against the WCS it should have recovered rather than
against another solver's opinion.

Two details make this worth having alongside the real-data tests in
``test_real_data.py``. It runs with no download and no network. And the synthetic WCS
(``conftest.make_wcs``) uses the *conventional* north-up/east-left parity, whereas a
real Seestar frame is mirrored -- so between them the two modules cover both entries of
``astrometry.SEED_PARITY``, and neither can quietly become the only one that works.

The star grid is deliberately jittered here. ``conftest.truth_table`` lays stars on an
exact 48-pixel lattice, which is close to the worst possible input for triangle matching:
every asterism has hundreds of congruent twins.
"""

import numpy as np
import pytest

pytest.importorskip("astroalign")

import conftest  # noqa: E402

from seestar_photometry import astrometry, frames  # noqa: E402


def jittered_truth(seed=17, jitter=11.0):
    """The standard truth table with positions perturbed off the lattice."""
    truth = conftest.truth_table()
    rng = np.random.default_rng(seed)
    truth["x"] = np.asarray(truth["x"]) + rng.uniform(-jitter, jitter, len(truth))
    truth["y"] = np.asarray(truth["y"]) + rng.uniform(-jitter, jitter, len(truth))
    return truth


@pytest.fixture(scope="module")
def truth():
    return jittered_truth()


@pytest.fixture(scope="module")
def catalogue(truth):
    return conftest.make_catalogue(truth)


def make_frame(directory, truth, name="synth.fit", solved=False, header=None):
    """A frame on disk with no ``.wcs`` sidecar, so a solve has real work to do."""
    path = conftest.write_cube(
        directory / name, conftest.make_cube(truth),
        header=header if header is not None else conftest.make_header("cube",
                                                                      solved=solved),
    )
    return frames.load_frame(path)


def wcs_disagreement(a, b, shape=(conftest.NY, conftest.NX)):
    """Worst positional disagreement between two WCS over the frame, in arcsec."""
    ny, nx = shape
    yy, xx = np.mgrid[0:ny:32, 0:nx:32]
    from seestar_photometry.gaiadb import separation_deg

    ra_a, dec_a = a.all_pix2world(xx.ravel(), yy.ravel(), 0)
    ra_b, dec_b = b.all_pix2world(xx.ravel(), yy.ravel(), 0)
    return max(separation_deg(p, q, r, s) * 3600.0
               for p, q, r, s in zip(ra_a, dec_a, ra_b, dec_b))


# --- the two routes to a solution ----------------------------------------------------

def test_solves_from_the_pointing_alone(tmp_path, truth, catalogue):
    """No header WCS at all -- the asterism bootstrap has to carry it.

    This is the case for anything ``stacking.stack_frame`` produced, and for every
    native Seestar stack, which carry a pointing but no WCS.
    """
    frame = make_frame(tmp_path, truth, "bootstrap.fit", solved=False)
    assert astrometry._header_wcs(frame) is None, "fixture must have no header WCS"

    wcs = astrometry.solve_local(frame, catalogue)
    assert wcs.has_celestial
    assert wcs_disagreement(wcs, conftest.make_wcs()) < conftest.PIXSCALE / 2


def test_solves_from_a_header_wcs(tmp_path, truth, catalogue):
    """The fast path. Cheaper, and the common case for a plate-solved exporter."""
    frame = make_frame(tmp_path, truth, "fastpath.fit", solved=True)
    assert astrometry._header_wcs(frame) is not None

    wcs = astrometry.solve_local(frame, catalogue)
    assert wcs_disagreement(wcs, conftest.make_wcs()) < conftest.PIXSCALE / 2


def test_a_wrong_header_wcs_falls_through_to_the_bootstrap(tmp_path, truth, catalogue):
    """A header solution that is wrong by degrees must not poison the result.

    The fast path is a guess, not a trusted input: if too few sources pair up it has to
    be abandoned rather than refined into a confident wrong answer.
    """
    header = conftest.make_header("cube", solved=True)
    header["CRVAL1"] = float(header["CRVAL1"]) + 3.0
    header["CRVAL2"] = float(header["CRVAL2"]) - 2.0
    frame = make_frame(tmp_path, truth, "badheader.fit", header=header)

    wcs = astrometry.solve_local(frame, catalogue)
    assert wcs_disagreement(wcs, conftest.make_wcs()) < conftest.PIXSCALE / 2


def test_recovers_a_rotated_field(tmp_path, truth, catalogue):
    """Alt-Az frames rotate; a real bundled stack needed 33 degrees of it."""
    rotated = truth.copy()
    angle = np.radians(40.0)
    cx, cy = conftest.NX / 2.0, conftest.NY / 2.0
    dx, dy = np.asarray(truth["x"]) - cx, np.asarray(truth["y"]) - cy
    rotated["x"] = cx + dx * np.cos(angle) - dy * np.sin(angle)
    rotated["y"] = cy + dx * np.sin(angle) + dy * np.cos(angle)
    # Stars rotated out of the frame would be injected off-image; keep the ones inside.
    inside = ((rotated["x"] > 20) & (rotated["x"] < conftest.NX - 20)
              & (rotated["y"] > 20) & (rotated["y"] < conftest.NY - 20))
    rotated = rotated[inside]

    frame = make_frame(tmp_path, rotated, "rotated.fit")
    wcs = astrometry.solve_local(frame, conftest.make_catalogue(truth))

    # The catalogue is built from the unrotated truth, so the recovered WCS must undo
    # the rotation: check the stars land where the catalogue says they are.
    ra, dec = wcs.all_pix2world(np.asarray(rotated["x"]), np.asarray(rotated["y"]), 0)
    reference = conftest.make_wcs().all_pix2world(
        np.asarray(truth["x"])[inside], np.asarray(truth["y"])[inside], 0
    )
    from seestar_photometry.gaiadb import separation_deg

    sep = separation_deg(reference[0][0], reference[1][0], ra[:1], dec[:1]) * 3600.0
    assert sep[0] < conftest.PIXSCALE


# --- the sidecar ----------------------------------------------------------------------

def test_caches_a_sidecar_and_reuses_it(tmp_path, truth, catalogue):
    frame = make_frame(tmp_path, truth, "cached.fit")
    cache = astrometry.wcs_cache_path(frame.path)
    assert not cache.exists()

    astrometry.solve_local(frame, catalogue)
    assert cache.exists()
    stamp = cache.stat().st_mtime_ns

    astrometry.solve_local(frame, catalogue)
    assert cache.stat().st_mtime_ns == stamp, "a cached solve must not rewrite"

    astrometry.solve_local(frame, catalogue, force=True)
    assert astrometry.load_wcs(frame).has_celestial


def test_reached_through_the_solve_dispatch(tmp_path, truth, catalogue):
    frame = make_frame(tmp_path, truth, "dispatch.fit")
    wcs = astrometry.solve(frame, solver="local", catalogue=catalogue)
    assert wcs_disagreement(wcs, conftest.make_wcs()) < conftest.PIXSCALE / 2


# --- failure is loud -------------------------------------------------------------------

def test_no_catalogue_says_how_to_supply_one(tmp_path, truth):
    frame = make_frame(tmp_path, truth, "nocat.fit")
    with pytest.raises(ValueError, match="catalogue"):
        astrometry.solve(frame, solver="local")


def test_unknown_solver_lists_local(tmp_path, truth):
    frame = make_frame(tmp_path, truth, "unknown.fit")
    with pytest.raises(ValueError, match="local"):
        astrometry.solve(frame, solver="nope")
    with pytest.raises(ValueError, match="local"):
        astrometry.solve_from_sources(frame, [1.0], [1.0], solver="nope")


def test_a_catalogue_for_the_wrong_field_raises(tmp_path, truth, catalogue):
    """Rather than caching a confident, wrong WCS."""
    header = conftest.make_header("cube")
    header["RA"] = float(header["RA"]) + 40.0
    frame = make_frame(tmp_path, truth, "wrongfield.fit", header=header)

    with pytest.raises(RuntimeError, match="could not match"):
        astrometry.solve_local(frame, catalogue)
    assert not astrometry.wcs_cache_path(frame.path).exists()


def test_a_frame_with_no_pointing_names_the_alternative(tmp_path, truth, catalogue):
    header = conftest.make_header("cube")
    del header["RA"]
    del header["DEC"]
    frame = make_frame(tmp_path, truth, "nopointing.fit", header=header)

    with pytest.raises(RuntimeError, match="astap"):
        astrometry.solve_local(frame, catalogue)


def test_too_few_sources_does_not_produce_a_wcs(tmp_path, truth, catalogue):
    """A three-star frame cannot be solved, and must say so rather than guess."""
    frame = make_frame(tmp_path, truth[:3], "sparse.fit")
    with pytest.raises(RuntimeError):
        astrometry.solve_local(frame, catalogue)
    assert not astrometry.wcs_cache_path(frame.path).exists()
