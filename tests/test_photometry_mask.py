"""Extended emission must not size the aperture.

The failure these guard against was found on a real M27 field: of 5469 green detections,
578 were bright and round, and exactly **one** cleared the 40 px isolation cut -- M27
itself (a = 42 px, b/a = 0.82). Isolation selects *for* extended objects in a crowded
field, because only something large has no close neighbour, and roundness does not exclude
a nebula. The curve of growth was then measured on the nebula and returned a 19.0 px
aperture where the stars wanted 5.3 px, with nothing reporting anything wrong.
"""

import numpy as np
import pytest

from seestar_photometry import frames, photometry

from conftest import FWHM, NOISE, NX, NY, make_cube, make_header, make_wcs


def _blob(cube, cx, cy, sigma, peak):
    """Add a big round Gaussian -- a stand-in for a planetary nebula."""
    yy, xx = np.mgrid[0:NY, 0:NX]
    g = peak * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sigma ** 2))
    for b in range(cube.shape[0]):
        cube[b] += g
    return cube


@pytest.fixture
def nebulous_frame(tmp_path):
    """A normal star field with one large bright blob dropped into a corner."""
    cube = make_cube()
    # Far from the star grid (MARGIN=42, PITCH=48) so it cannot blend with a real star,
    # and wide enough that its semi-major axis is many times the stellar one.
    _blob(cube, cx=NX - 40, cy=NY - 40, sigma=9.0, peak=40.0 * NOISE)
    return frames.SeestarFrame(
        data=cube, header=make_header(), model="S50",
        path=tmp_path / "nebulous.fits", layout="cube",
    )


def test_nebula_is_rejected_from_the_curve_of_growth(nebulous_frame):
    """The size cut keeps the blob out, so the aperture stays stellar."""
    cog = photometry.curve_of_growth(nebulous_frame, band="G")
    assert cog.meta["n_stars"] >= photometry.MIN_COG_STARS
    ap = photometry.aperture_for_enclosed_flux(cog, 0.90)
    # The Gaussian expectation, matching test_aperture_radius_matches_gaussian_expectation.
    assert ap == pytest.approx(0.911 * FWHM["G"], rel=0.15)


def test_size_cut_is_what_rejects_it(nebulous_frame):
    """Disabling the size cut lets the blob into the sample -- the regression itself."""
    strict = photometry.curve_of_growth(nebulous_frame, band="G")
    loose = photometry.curve_of_growth(nebulous_frame, band="G", size_ratio=np.inf)
    assert loose.meta["n_stars"] > strict.meta["n_stars"]
    # and it drags the aperture outwards
    assert (photometry.aperture_for_enclosed_flux(loose, 0.90)
            > photometry.aperture_for_enclosed_flux(strict, 0.90))


def test_mask_excludes_a_region_from_detection(nebulous_frame):
    """Detections whose centroid lands on a masked pixel are dropped."""
    mask = np.zeros((NY, NX), dtype=bool)
    mask[NY - 90:, NX - 90:] = True
    ext_all = photometry.extract_sources(nebulous_frame, n_fwhm=2.0)
    ext_masked = photometry.extract_sources(nebulous_frame, n_fwhm=2.0, mask=mask)
    x = np.asarray(ext_masked.band("G")["x"], float)
    y = np.asarray(ext_masked.band("G")["y"], float)
    assert not np.any((x > NX - 90) & (y > NY - 90))
    assert len(ext_masked.band("G")) < len(ext_all.band("G"))


def test_mask_does_not_change_unmasked_fluxes(nebulous_frame):
    """Masking must not silently alter the flux of a star far from the region."""
    mask = np.zeros((NY, NX), dtype=bool)
    mask[NY - 90:, NX - 90:] = True
    a = photometry.extract_sources(nebulous_frame, aperture=5.0)
    b = photometry.extract_sources(nebulous_frame, aperture=5.0, mask=mask)
    ga, gb = a.band("G"), b.band("G")
    # match on position, then compare flux for the stars present in both
    ka = {(round(float(r["x"])), round(float(r["y"]))): float(r["flux"]) for r in ga}
    kb = {(round(float(r["x"])), round(float(r["y"]))): float(r["flux"]) for r in gb}
    shared = [k for k in ka if k in kb and k[0] < NX - 120 and k[1] < NY - 120]
    assert len(shared) > 10
    for k in shared:
        # background changes very slightly (the blob no longer biases it), so allow 1%
        assert kb[k] == pytest.approx(ka[k], rel=0.01)


def test_few_cog_stars_falls_back_to_fwhm(nebulous_frame):
    """An unusable COG must yield the 1.2 x FWHM fallback, not a one-star radius."""
    # An impossible isolation requirement empties the sample.
    ext = photometry.extract_sources(nebulous_frame, enclosed=0.90, isolation=10_000)
    fwhm = ext.fwhm[1]
    assert ext.aperture[1] == pytest.approx(1.2 * fwhm, rel=1e-6)


def test_empty_cog_is_reported_not_warned(nebulous_frame):
    """No usable stars -> all-nan curve, n_stars=0, and no numpy warning."""
    with np.errstate(all="raise"):
        cog = photometry.curve_of_growth(nebulous_frame, band="G", isolation=10_000)
    assert cog.meta["n_stars"] == 0
    assert np.all(~np.isfinite(np.asarray(cog["flux_frac"], float)))


def test_sky_mask_places_a_circle_from_sky_coordinates():
    wcs = make_wcs()
    ra, dec = wcs.pixel_to_world_values(NX / 2, NY / 2)
    mask = photometry.sky_mask((NY, NX), wcs, float(ra), float(dec), 60.0)
    assert mask.shape == (NY, NX)
    assert mask[int(NY / 2), int(NX / 2)]
    assert not mask[0, 0]
    # a zero-radius region masks nothing
    assert not photometry.sky_mask((NY, NX), wcs, float(ra), float(dec), 0.0).any()


def test_mask_shape_is_validated(nebulous_frame):
    with pytest.raises(ValueError, match="mask shape"):
        photometry.extract_sources(nebulous_frame, mask=np.zeros((5, 5), bool))
