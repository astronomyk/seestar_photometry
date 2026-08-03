"""Reference-catalogue handling: a cached oversized Gaia field, and cross-match.

A Seestar dithers through a night, so a dataset covers more sky than any single
FOV. We therefore build one **oversized catalogue once** -- tiling the area into
overlapping cones queried in parallel, then concatenating and de-duplicating -- and
cache it as an ECSV. Every frame then subsets the cached catalogue to its own
footprint (:func:`sources_in_frame`); no per-frame catalogue queries are needed,
which is what makes the per-frame work pure CPU and trivially parallel.

The query carries Gaia DR3 synthetic Johnson-Kron-Cousins V (``v_jkc_mag``) -- the
reference magnitude for this project -- plus ``b_jkc_mag``/``r_jkc_mag`` for a
synthetic B-R colour (left join on ``gaiadr3.synthetic_photometry_gspc``, null where
absent), and ``phot_variable_flag`` for vetting variables later. Gaia synthetic V
beats APASS here: it is homogeneous, and its colour information is what makes the
per-frame colour term fittable.

See ``docs/astrometry-and-gaia.md``.
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.table import MaskedColumn, Table, unique, vstack

#: Catalogue columns carried through the field query and the cross-match.
GAIA_COLUMNS = (
    "source_id", "ra", "dec",
    "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp",
    "phot_variable_flag",
    "v_jkc_mag", "b_jkc_mag", "r_jkc_mag",
)


def _cone_query(ra, dec, radius, gmag_limit):
    return f"""
        SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag,
               g.phot_bp_mean_mag, g.phot_rp_mean_mag, g.bp_rp,
               g.phot_variable_flag,
               s.v_jkc_mag, s.b_jkc_mag, s.r_jkc_mag
        FROM gaiadr3.gaia_source AS g
        LEFT OUTER JOIN gaiadr3.synthetic_photometry_gspc AS s
            ON g.source_id = s.source_id
        WHERE 1 = CONTAINS(POINT('ICRS', g.ra, g.dec),
                           CIRCLE('ICRS', {ra}, {dec}, {radius}))
          AND g.phot_g_mean_mag < {gmag_limit}
    """


def _fetch_cone(ra, dec, radius, gmag_limit):
    """Run a single Gaia cone search. Fresh client per call for thread safety."""
    from astroquery.gaia import GaiaClass

    gaia = GaiaClass()
    gaia.ROW_LIMIT = -1
    return gaia.launch_job_async(_cone_query(ra, dec, radius, gmag_limit)).get_results()


def _tile_centers(center, half_size_deg, n_tiles):
    """Grid of tile centres (ra, dec) as true angular offsets around ``center``.

    Offsets are taken in a :class:`SkyOffsetFrame` rather than by adding degrees to
    RA, so the tiling stays correct near the pole -- which matters, since MW Cam sits
    at Dec +81 where naive RA steps collapse.
    """
    origin = SkyCoord(center[0] * u.deg, center[1] * u.deg)
    frame = SkyOffsetFrame(origin=origin)
    cell = 2 * half_size_deg / n_tiles
    offsets = [-half_size_deg + cell * (i + 0.5) for i in range(n_tiles)]
    centers = []
    for dx in offsets:
        for dy in offsets:
            p = SkyCoord(lon=dx * u.deg, lat=dy * u.deg, frame=frame).icrs
            centers.append((float(p.ra.deg), float(p.dec.deg)))
    return centers, cell


def fetch_gaia_mosaic(
    center, cache_path, half_size_deg=1.5, n_tiles=1, gmag_limit=17.0,
    max_workers=5, overwrite=False, attempts=4,
):
    """Build (and cache) an oversized Gaia catalogue by tiling a box.

    Covers a ``2 * half_size_deg`` square centred on ``center`` with an
    ``n_tiles x n_tiles`` grid of overlapping cones, then concatenates and
    de-duplicates on ``source_id``. Cached as an ECSV and read back on later calls.

    The Gaia TAP endpoint intermittently truncates an async result
    (``IncompleteRead``), so the whole build is retried a few times. The cache is
    only written on full success -- a partial mosaic silently missing a tile is much
    worse than a failed build, because every downstream frame in that region would
    lose its calibration stars without any error.

    Parameters
    ----------
    center : tuple of float
        Field centre ``(ra, dec)`` in degrees.
    cache_path : path-like
        ECSV cache for the full catalogue.
    half_size_deg : float
        Half-size of the box. 1.5 deg comfortably covers a night's dithering.
    n_tiles : int
        Tiles per axis. Keep at 1 unless the box is large: Gaia TAP is unreliable
        under concurrent jobs, and a single cone is both faster and safer.
    """
    cache_path = Path(cache_path)
    if cache_path.exists() and not overwrite:
        return Table.read(cache_path)

    import time

    centers, cell = _tile_centers(center, half_size_deg, n_tiles)
    radius = cell / 2 * np.sqrt(2) + 0.1  # cover each cell (incl. corners) + overlap
    last = None
    for attempt in range(attempts):
        try:
            with ThreadPoolExecutor(max_workers=max_workers) as ex:
                tables = list(
                    ex.map(lambda c: _fetch_cone(c[0], c[1], radius, gmag_limit), centers)
                )
            catalogue = unique(vstack(tables, metadata_conflicts="silent"), keys="source_id")
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            catalogue.write(cache_path, overwrite=True)
            return catalogue
        except Exception as exc:  # transient TAP truncation / dropped connection
            last = exc
            if attempt < attempts - 1:
                time.sleep(5 * (attempt + 1))
    raise last


def load_catalogue(cache_path):
    """Read a cached reference-catalogue ECSV."""
    return Table.read(cache_path)


def sources_in_frame(catalogue, wcs, frame, margin=0.0):
    """Subset a catalogue to the sources falling within a frame.

    Maps catalogue sky positions to pixel coordinates via ``wcs`` and keeps those
    inside the image (optionally expanded by ``margin`` pixels). The shape comes
    from the frame data, so both FITS layouts work.
    """
    ny, nx = frame.shape
    x, y = wcs.world_to_pixel_values(
        np.asarray(catalogue["ra"]), np.asarray(catalogue["dec"])
    )
    inside = (x >= -margin) & (x < nx + margin) & (y >= -margin) & (y < ny + margin)
    return catalogue[inside]


def crossmatch_table(sources, wcs, catalogue, tol_arcsec=2.0, columns=GAIA_COLUMNS):
    """Augment a sources table with sky position and matched catalogue columns.

    Every input row is kept and gains ``ra``, ``dec`` (from the WCS) and
    ``sep_arcsec`` (distance to its nearest catalogue source). The catalogue columns
    are filled where the match is within ``tol_arcsec`` and **masked** otherwise --
    including where the matched source itself lacks the value (e.g. no synthetic V),
    so a non-masked ``v_jkc_mag`` always means a usable match. Returns a new table.
    """
    ra, dec = wcs.all_pix2world(np.asarray(sources["x"]), np.asarray(sources["y"]), 0)
    det = SkyCoord(ra * u.deg, dec * u.deg)
    cat = SkyCoord(
        np.asarray(catalogue["ra"]) * u.deg, np.asarray(catalogue["dec"]) * u.deg
    )
    idx, d2d, _ = det.match_to_catalog_sky(cat)
    unmatched = d2d.arcsec >= tol_arcsec

    out = sources.copy()
    out["ra"] = ra
    out["dec"] = dec
    out["sep_arcsec"] = d2d.arcsec
    for col in columns:
        if col not in catalogue.colnames:
            continue
        name = {"ra": "cat_ra", "dec": "cat_dec"}.get(col, col)
        values = np.ma.getdata(catalogue[col])[idx]
        mask = unmatched | np.ma.getmaskarray(catalogue[col])[idx]
        out[name] = MaskedColumn(values, mask=mask)
    return out
