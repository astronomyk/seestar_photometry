#!/usr/bin/env python
"""Build a local Gaia catalogue in the :mod:`gaiadb` layout.

Two scales, and only the first is implemented here:

**A region**, queried straight from Gaia TAP::

    uv run python tools/build_gaia_catalogue.py --ra 186.6821 --dec 81.474 --radius 3

This is what you want for trying the offline path on a field you already have data for,
and for the real-data tests. A few degrees is a comfortable query; a few tens of degrees
is not, because TAP will truncate it.

**The whole sky** is a different job and is deliberately not attempted through TAP. At
G < 17.65 the synthetic-photometry catalogue is ~220 million rows, which is a bulk
download from the Gaia archive (the per-HEALPix ``gaia_source`` and
``synthetic_photometry_gspc`` files) joined offline, not an ADQL query. The writing half
is already shared -- :func:`gaiadb.write_dataset` takes whatever table it is handed --
so the remaining work is the fetch-and-join, not the layout.

Maintenance script, not shipped in the wheel.

Output goes to ``--out`` (default: a directory named for the dataset under the cache
:mod:`gaiadb` reads from, so a build is immediately live).
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from seestar_photometry import gaiadb  # noqa: E402

#: Query template. An INNER join, not a left outer one: a source with no synthetic
#: photometry cannot be a calibrator and cannot be cut on V, so it has no place in a
#: dataset whose whole purpose is calibrators.
#:
#: ``c_star`` and ``v_jkc_flag`` come from the GSPC table -- ``c_star`` is the corrected
#: BP/RP flux-excess statistic (the blend indicator) and ``v_jkc_flag`` says whether the
#: source sits inside the range the synthetic V was validated over.
QUERY = """
    SELECT g.source_id, g.ra, g.dec,
           g.phot_g_mean_mag, g.bp_rp, g.phot_variable_flag,
           s.v_jkc_mag, s.b_jkc_mag, s.r_jkc_mag, s.v_jkc_flag,
           g.pmra, g.pmdec,
           s.c_star, g.ruwe, g.teff_gspphot,
           g.ipd_frac_multi_peak, g.non_single_star,
           g.duplicated_source, g.in_galaxy_candidates, g.has_epoch_photometry
    FROM gaiadr3.gaia_source AS g
    JOIN gaiadr3.synthetic_photometry_gspc AS s ON g.source_id = s.source_id
    WHERE 1 = CONTAINS(POINT('ICRS', g.ra, g.dec),
                       CIRCLE('ICRS', {ra}, {dec}, {radius}))
      AND s.v_jkc_mag < {v_limit}
"""


def fetch(ra, dec, radius, v_limit):
    """One TAP cone, returning the raw astropy Table."""
    from astroquery.gaia import GaiaClass

    gaia = GaiaClass()
    gaia.ROW_LIMIT = -1
    query = QUERY.format(ra=ra, dec=dec, radius=radius, v_limit=v_limit)
    print(f"[gaia] querying {radius} deg around ({ra}, {dec}), V < {v_limit}",
          flush=True)
    return gaia.launch_job_async(query).get_results()


def conform(table):
    """Coerce a TAP result into the :mod:`gaiadb` schema.

    Mostly a dtype pass. The one judgement call is the boolean columns: TAP hands them
    back as objects or as ``'t'``/``'f'`` strings depending on the client version, and a
    naive ``astype(bool)`` on the string form makes *everything* True.
    """
    from astropy.table import MaskedColumn, Table

    out = Table()
    for name in gaiadb.COLUMNS:
        if name not in table.colnames:
            raise SystemExit(
                f"the TAP result has no {name!r} column. The archive's schema may have "
                f"changed; check the datamodel and update QUERY.\n  got: "
                f"{sorted(table.colnames)}"
            )
        column = table[name]
        mask = np.ma.getmaskarray(column)
        values = np.ma.getdata(column)
        dtype = gaiadb.DTYPES[name]

        if dtype == "bool":
            values = _as_bool(values)
        elif dtype == "str":
            values = np.asarray(values).astype(str)
        elif np.issubdtype(np.dtype(dtype), np.integer):
            values = np.nan_to_num(np.asarray(values, dtype="float64"), nan=0.0)
            values = values.astype(dtype)
        else:
            values = np.asarray(values, dtype=dtype)
        out[name] = MaskedColumn(values, mask=mask, name=name,
                                 unit=gaiadb.UNITS.get(name))
    return out


def _as_bool(values):
    """Booleans out of whatever TAP produced -- real bools, or 't'/'f', or 0/1."""
    values = np.asarray(values)
    if values.dtype == bool:
        return values
    if values.dtype.kind in "OSU":
        text = np.char.lower(values.astype(str))
        return np.isin(text, ("true", "t", "1", "yes"))
    return values.astype(bool)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ra", type=float, required=True, help="cone centre RA, degrees")
    ap.add_argument("--dec", type=float, required=True, help="cone centre Dec, degrees")
    ap.add_argument("--radius", type=float, default=3.0, help="cone radius, degrees")
    ap.add_argument("--v-limit", type=float, default=gaiadb.V_LIMIT,
                    help="faint limit in synthetic V")
    ap.add_argument("--out", type=Path, default=None,
                    help="output directory (default: the gaiadb cache, so it goes live)")
    args = ap.parse_args()

    if args.radius > 10:
        print(f"[gaia] warning: {args.radius} deg is a large single cone; TAP truncates "
              "these. Build several smaller ones into the same --out instead.",
              flush=True)

    raw = fetch(args.ra, args.dec, args.radius, args.v_limit)
    print(f"[gaia] {len(raw)} rows", flush=True)
    if not len(raw):
        raise SystemExit("no rows returned; check the coordinates")

    table = conform(raw)
    out = Path(args.out) if args.out else gaiadb.dataset_dir()

    # Merge rather than replace: building several cones into one directory is how a
    # multi-field or larger-than-one-query dataset gets made. Parts are whole files, so
    # a pixel touched by two cones has to be rebuilt from the union of both.
    existing = _existing_rows(out, table)
    if existing is not None:
        from astropy.table import unique, vstack

        before = len(table)
        table = unique(vstack([existing, table], metadata_conflicts="silent"),
                       keys="source_id")
        print(f"[gaia] merged with {len(existing)} rows already present "
              f"({before} new, {len(table)} total)", flush=True)

    meta = gaiadb.write_dataset(table, out)
    total = sum(p["bytes"] for p in meta["parts"].values())
    print(f"[gaia] wrote {meta['rows']} rows into {len(meta['parts'])} parts, "
          f"{total / 1e6:.1f} MB ({total / max(meta['rows'], 1):.1f} B/row)")
    print(f"[gaia] -> {out}")
    if args.out is None:
        print("[gaia] this is the directory gaiadb reads, so it is live now; "
              "Project(catalogue_backend='auto') will pick it up")


def _existing_rows(directory, incoming):
    """Rows already stored in the parts the new table would overwrite, or ``None``."""
    directory = Path(directory)
    if gaiadb.manifest(directory) is None:
        return None
    pixels = np.unique(gaiadb.hpx_of(incoming["source_id"]))
    paths = [p for p in (gaiadb.part_path(i, directory) for i in pixels) if p.exists()]
    if not paths:
        return None

    import pyarrow as pa
    import pyarrow.parquet as pq

    arrow = pa.concat_tables([pq.read_table(p, columns=list(gaiadb.COLUMNS))
                              for p in paths])
    return gaiadb._to_table(arrow)


if __name__ == "__main__":
    main()
