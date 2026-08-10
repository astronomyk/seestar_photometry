#!/usr/bin/env python
"""Build the all-sky Gaia catalogue by joining GSPC photometry onto existing positions.

The obvious route -- bulk-download everything -- does not work here, and the reason is
worth writing down because it is counter-intuitive.

What the catalogue needs comes from two Gaia tables. ``synthetic_photometry_gspc`` has
the synthetic Johnson V and the colour, and bulk-downloads as 44 files totalling ~42 GB.
``gaia_source`` has the positions and proper motions, and bulk-downloads as 3387 files
totalling **~790 GB** -- 152 columns for 1.8 billion rows, to keep 6 columns for the 220
million that have synthetic photometry. Neither table alone is enough: GSPC carries no
``ra``/``dec`` at all, only ``source_id``.

Nor can TAP do the join in ranges. Chunking on ``source_id BETWEEN`` has no spatial
index behind it, and the archive answers with an HTTP 500 after exactly 182 seconds --
measured, three times, on queries as simple as a row count.

So this joins **new photometry onto positions already in hand**. Given a directory of
per-HEALPix parquet tiles carrying ``source_id``/``ra``/``dec``/``pmra``/``pmdec``, only
the 42 GB of GSPC has to be fetched.

Two passes:

1. **Extract.** Stream each GSPC file, keep five columns, and append rows to bucket
   files by HEALPix. The GSPC files are not sorted by ``source_id`` -- rows scatter
   across the whole sky within one file -- so a merge join is not available and the
   rows have to be routed by destination first.
2. **Join.** For each bucket, index it by ``source_id``, join onto that range of tiles,
   apply the V cut, and write the :mod:`gaiadb` dataset.

Both passes are resumable: a completed GSPC file leaves a marker, and pass 2 skips parts
already written. Kill it and run it again.

Usage::

    uv run python tools/build_gaia_allsky.py \\
        --tiles D:/tmp/gaia_V/nside32 --work D:/tmp/gaia_build --out D:/tmp/gaia_allsky

Maintenance script, not shipped in the wheel. See ``tools/HOSTING.md``.
"""

import argparse
import shutil
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from seestar_photometry import gaiadb  # noqa: E402

#: Where the GSPC bulk files live, and how many there are.
GSPC_URL = ("http://cdn.gea.esac.esa.int/Gaia/gdr3/Performance_verification/"
            "synthetic_photometry_gspc")
GSPC_FILES = 44

#: The only GSPC columns worth carrying. ``v_jkc_mag`` is fetched too, as a cross-check
#: against the value already in the tiles rather than because it is needed.
GSPC_COLUMNS = ["source_id", "v_jkc_mag", "b_jkc_mag", "r_jkc_mag", "v_jkc_flag",
                "c_star"]

#: Columns expected in the existing position tiles.
TILE_COLUMNS = ["source_id", "ra", "dec", "pmra", "pmdec", "phot_g_mean_mag",
                "v_jkc_mag"]

#: HEALPix tiles per bucket. 48 buckets over 12288 tiles keeps each one to ~100 MB,
#: which is the point: the join is done a bucket at a time so memory stays flat.
TILES_PER_BUCKET = 256

#: Fixed-width record written to the bucket files. Raw numpy rather than parquet
#: because pass 1 appends to 48 files interleaved and this is the cheapest thing that
#: can be appended to.
BUCKET_DTYPE = np.dtype([
    ("source_id", "<i8"), ("v_jkc_mag", "<f4"), ("b_jkc_mag", "<f4"),
    ("r_jkc_mag", "<f4"), ("c_star", "<f4"), ("v_jkc_flag", "u1"),
])


def bucket_of(hpx):
    return hpx // TILES_PER_BUCKET


def n_buckets():
    return -(-gaiadb.N_PARTS // TILES_PER_BUCKET)


# --- pass 1: fetch and route -------------------------------------------------------------

def _header_offset(path):
    """Number of leading ``#`` comment lines before the CSV header.

    Every Gaia bulk file starts with a block of provenance comments, and pyarrow's CSV
    reader has no notion of a comment character -- it would take the first ``#`` line as
    the header and then fail to find any column by name. Counted per file rather than
    hardcoded, because the block's length is not part of any contract.
    """
    import gzip

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for n, line in enumerate(fh):
            if not line.startswith("#"):
                return n
    raise ValueError(f"{path} has no header line")


def _marker(index, work):
    return work / "markers" / f"gspc_{index}.done"


def download_one(index, work, quiet=False):
    """Fetch one GSPC file, returning its local path. Safe to run concurrently."""
    name = f"SyntheticPhotometryGspc_{index}.csv.gz"
    local = work / "gspc" / name
    local.parent.mkdir(parents=True, exist_ok=True)
    if local.exists():
        return local
    t0 = time.time()
    tmp = local.with_suffix(f".{index}.part")
    with urllib.request.urlopen(f"{GSPC_URL}/{name}", timeout=600) as response, \
            open(tmp, "wb") as fh:
        shutil.copyfileobj(response, fh, length=1 << 22)
    tmp.replace(local)
    if not quiet:
        mb = local.stat().st_size / 1e6
        print(f"  [{index:2d}] fetched {mb:6.0f} MB in {time.time()-t0:5.0f}s "
              f"({mb/max(time.time()-t0, 1):.1f} MB/s)", flush=True)
    return local


def extract_all(work, workers=4, keep_downloads=False, quiet=False):
    """Fetch every GSPC file and route its rows into the buckets.

    Downloads run concurrently and routing does not, which is the shape the work wants:
    the CDN gives ~2 MB/s on one stream and ~5 MB/s over four, while routing a whole
    file takes 16 seconds. So the download is the wall clock and the four-way pool
    roughly halves it, while bucket appends stay single-threaded and need no locking.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    todo = [i for i in range(GSPC_FILES) if not _marker(i, work).exists()]
    if not quiet:
        print(f"  {len(todo)} of {GSPC_FILES} files to fetch "
              f"({GSPC_FILES - len(todo)} already routed)", flush=True)
    total = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_one, i, work, quiet): i for i in todo}
        for future in as_completed(futures):
            index = futures[future]
            total += route_one(index, future.result(), work, keep_downloads, quiet)
    return total


def route_one(index, local, work, keep_download=False, quiet=False):
    """Append one downloaded GSPC file's rows to the bucket files."""
    import pyarrow.csv as pacsv

    marker = _marker(index, work)
    if marker.exists():
        return 0

    t0 = time.time()
    table = pacsv.read_csv(
        local,
        read_options=pacsv.ReadOptions(block_size=1 << 26,
                                       skip_rows=_header_offset(local)),
        parse_options=pacsv.ParseOptions(),
        convert_options=pacsv.ConvertOptions(include_columns=GSPC_COLUMNS),
    )
    record = np.empty(table.num_rows, dtype=BUCKET_DTYPE)
    record["source_id"] = np.asarray(table["source_id"])
    for column in ("v_jkc_mag", "b_jkc_mag", "r_jkc_mag", "c_star"):
        # Nulls arrive as NaN, which is what the join wants: a missing magnitude has to
        # stay missing rather than become a zero that reads as an extremely bright star.
        record[column] = np.asarray(table[column], dtype="f4")
    record["v_jkc_flag"] = np.nan_to_num(
        np.asarray(table["v_jkc_flag"], dtype="f4"), nan=0.0
    ).astype("u1")

    hpx = record["source_id"] >> 49
    buckets = bucket_of(hpx)
    order = np.argsort(buckets, kind="stable")
    record, buckets = record[order], buckets[order]
    edges = np.searchsorted(buckets, np.arange(n_buckets() + 1))
    written = 0
    for b in range(n_buckets()):
        chunk = record[edges[b]:edges[b + 1]]
        if not len(chunk):
            continue
        path = work / "buckets" / f"b{b:04d}.bin"
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "ab") as fh:
            chunk.tofile(fh)
        written += len(chunk)

    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(f"{written}\n", encoding="utf-8")
    if not keep_download:
        local.unlink(missing_ok=True)
    if not quiet:
        print(f"  [{index:2d}] routed {written/1e6:5.1f}M rows in "
              f"{time.time()-t0:5.0f}s", flush=True)
    return written


# --- pass 2: join ---------------------------------------------------------------------

def join_bucket(bucket, tiles_dir, work, out, v_limit, quiet=False):
    """Join one bucket's photometry onto its tiles and write those parts."""
    import pyarrow.parquet as pq
    from astropy.table import MaskedColumn, Table

    path = work / "buckets" / f"b{bucket:04d}.bin"
    if not path.exists():
        return 0, 0
    photometry = np.fromfile(path, dtype=BUCKET_DTYPE)
    order = np.argsort(photometry["source_id"], kind="stable")
    photometry = photometry[order]
    ids = photometry["source_id"]

    rows = written = 0
    for tile in range(bucket * TILES_PER_BUCKET,
                      min((bucket + 1) * TILES_PER_BUCKET, gaiadb.N_PARTS)):
        source = tiles_dir / f"tile_{tile:05d}.parquet"
        if not source.exists():
            continue
        part = gaiadb.part_path(tile, out)
        if part.exists():
            continue

        have = pq.read_table(source, columns=TILE_COLUMNS)
        sid = have["source_id"].to_numpy()
        if not len(sid):
            continue

        # Positions are the spine; photometry is looked up per source. A source with no
        # GSPC row cannot be a calibrator and is dropped by the V cut below anyway.
        where = np.searchsorted(ids, sid)
        where = np.clip(where, 0, len(ids) - 1)
        hit = ids[where] == sid

        v = np.where(hit, photometry["v_jkc_mag"][where], np.nan)
        keep = np.isfinite(v) & (v < v_limit)
        if not keep.any():
            continue
        idx = where[keep]

        table = Table()
        table["source_id"] = sid[keep]
        for column in ("ra", "dec"):
            table[column] = have[column].to_numpy()[keep]
        for column in ("pmra", "pmdec", "phot_g_mean_mag"):
            values = np.asarray(have[column].to_numpy(zero_copy_only=False), dtype="f4")
            table[column] = MaskedColumn(np.nan_to_num(values[keep]),
                                         mask=~np.isfinite(values[keep]))
        for column in ("v_jkc_mag", "b_jkc_mag", "r_jkc_mag", "c_star"):
            values = photometry[column][idx]
            table[column] = MaskedColumn(np.nan_to_num(values),
                                         mask=~np.isfinite(values))
        table["v_jkc_flag"] = MaskedColumn(
            photometry["v_jkc_flag"][idx].astype(bool),
            mask=np.zeros(int(keep.sum()), dtype=bool),
        )
        gaiadb.write_part(table, out, pixel=tile)
        rows += len(table)
        written += 1
    if not quiet:
        print(f"  bucket {bucket:3d}: {written:4d} parts, {rows/1e6:5.2f}M rows",
              flush=True)
    return rows, written


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles", type=Path, required=True,
                    help="directory of tile_NNNNN.parquet position tiles")
    ap.add_argument("--work", type=Path, required=True, help="scratch for buckets")
    ap.add_argument("--out", type=Path, required=True, help="dataset to write")
    ap.add_argument("--v-limit", type=float, default=gaiadb.V_LIMIT)
    ap.add_argument("--keep-downloads", action="store_true",
                    help="do not delete each GSPC file after routing it (42 GB)")
    ap.add_argument("--workers", type=int, default=4,
                    help="concurrent downloads; the CDN gives ~2 MB/s per stream")
    ap.add_argument("--stage", choices=("extract", "join", "all"), default="all")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    if args.stage in ("extract", "all"):
        print(f"pass 1: {GSPC_FILES} GSPC files (~42 GB) -> {n_buckets()} buckets",
              flush=True)
        total = extract_all(args.work, workers=args.workers,
                            keep_downloads=args.keep_downloads)
        print(f"pass 1 done: {total/1e6:.1f}M rows routed in "
              f"{(time.time()-t_start)/60:.0f} min\n", flush=True)

    if args.stage in ("join", "all"):
        print(f"pass 2: joining onto {args.tiles}", flush=True)
        rows = parts = 0
        for b in range(n_buckets()):
            r, w = join_bucket(b, args.tiles, args.work, args.out, args.v_limit)
            rows += r
            parts += w
        gaiadb.finalise_manifest(args.out)
        print(f"\npass 2 done: {rows/1e6:.1f}M rows in {parts} parts")
        size = sum(p.stat().st_size for p in args.out.rglob("*.parquet"))
        print(f"dataset: {size/1e9:.2f} GB, {size/max(rows,1):.1f} B/row")
        print(f"total: {(time.time()-t_start)/60:.0f} min -> {args.out}")


if __name__ == "__main__":
    main()
