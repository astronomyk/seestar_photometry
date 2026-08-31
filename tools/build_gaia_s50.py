#!/usr/bin/env python
"""Add synthetic Seestar S50 magnitudes to the all-sky Gaia catalogue.

Integrates every Gaia DR3 XP sampled mean spectrum through the S50's fitted R/G/B
passbands and writes the six new columns (:data:`S50_COLUMNS`) into a copy of the
:mod:`gaiadb` dataset. Calibrating the Bayer planes against these instead of
``v_jkc_mag``/``b_jkc_mag``/``r_jkc_mag`` removes the zero point's colour term: measured
on 57 fifteen-minute 20 Cep stacks, G +0.054 -> -0.002, R -0.132 -> -0.002, B +0.548 ->
~0 mag/mag, and the B-plane calibrator rms fell from 347 to 104 mmag.

**The spectra reach G < 15.0 and no further.** ESA's bulk ``xp_sampled_mean_spectrum``
product is the bright subset of the XP population. Measured over the finished all-sky
build, coverage of the V < 17.5 catalogue is 100% below G = 15.0000 and 0% above it --
34,238,688 of 148,667,792 rows, 23.03%. So the six columns are populated for the bright
quarter of the catalogue and masked for the rest, which is the shape the calibration wants
anyway (a saturated-to-faint Seestar stack draws its calibrators from roughly V 11-16).
Full coverage would mean the ``xp_continuous_mean_spectrum`` product instead: 3.6 TB
against 114 GB, and a GaiaXPy sampling pass per source rather than a numpy dot product.

Data source: the **bulk CDN**, not the archive query services, which hung for 15+ minutes
on 100-source datalink requests on two separate evenings. 3386 gzipped ECSV files, one
per HEALPix-8 range, 114.4 GB total, tiling hp8 0-786431 contiguously.

The id-to-file mapping is the whole trick. A Gaia ``source_id`` encodes its nested
HEALPix index in the high bits, so ``source_id >> 43`` is its level-8 pixel just as
``source_id >> 49`` is the level-5 :mod:`gaiadb` partition. Each part is therefore
exactly 64 consecutive hp8 pixels and maps to a contiguous run of bulk files.

Two stages, both resumable:

1. **fetch.** Download the bulk files. Restarting skips whatever is already the right
   size on disk.
2. **synth.** Parallel over *groups of consecutive parts*. A group parses each of its
   files once, integrates, and writes its parts; restarting skips parts already written.
   Only the file straddling a group boundary is parsed twice, which is ~11% of a stage
   whose cost is dominated by the download anyway.

Usage::

    uv run python tools/build_gaia_s50.py --stage list --work D:/tmp/xp
    uv run python tools/build_gaia_s50.py --validate --work D:/tmp/xp
    uv run python tools/build_gaia_s50.py --src D:/tmp/gaia_allsky --out D:/tmp/gaia_s50 \\
        --work D:/tmp/xp --parts 3913-3919 --workers 6

Maintenance script, not shipped in the wheel. See ``tools/HOSTING.md``.
"""

import argparse
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.request
from pathlib import Path

# Before numpy, and therefore before anything can read them: the pool already uses every
# core, and a BLAS that starts its own threads per worker oversubscribes the machine.
# Setting them in the worker would be too late under fork, where numpy is already imported.
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from seestar_photometry import gaiadb  # noqa: E402

#: Where the bulk spectra are served from, and the S3-style listing behind the same
#: storage. The listing is paginated by ``marker``; four pages cover the product.
CDN = ("https://cdn.gea.esac.esa.int/Gaia/gdr3/Spectroscopy/"
       "xp_sampled_mean_spectrum/")
BUCKET = "https://gaia.eu-1.cdn77-storage.com/"
LISTING_PREFIX = ("Gaia/gdr3/Spectroscopy/xp_sampled_mean_spectrum/"
                  "XpSampledMeanSpectrum_")

#: Cached copy of the listing, so ``fetch`` and ``synth`` need no network agreement.
LISTING = "xp_listing.json"

#: Attempts per file before giving up on it. See :func:`fetch_one`.
FETCH_ATTEMPTS = 5

#: The product's fixed sampling: 336-1020 nm in 2 nm steps, 343 values.
XP_GRID = np.arange(336.0, 1021.0, 2.0)
N_WL = len(XP_GRID)

#: Non-spectrum columns ahead of ``flux`` in each row: source_id, solution_id, ra, dec.
N_LEAD = 4

C_NM = 2.99792458e17   # speed of light [nm/s]
F_AB = 3631e-26        # AB reference [W m^-2 Hz^-1]

#: HEALPix level the bulk files are partitioned by, and the shift that recovers it from
#: a ``source_id``. Gaia packs the level-12 nested index above 35 bits of running number.
HPX8_LEVEL = 8
HPX8_SHIFT = 35 + 2 * (12 - HPX8_LEVEL)

#: hp8 pixels inside one :mod:`gaiadb` part. 64, and the reason a part maps to a
#: contiguous run of files rather than a scatter of them.
HP8_PER_PART = 4 ** (HPX8_LEVEL - gaiadb.HPX_LEVEL)

#: The passband file, and the columns in it. Copied from ``photometry-adventures-2``
#: rather than rebuilt: it is the fitted instrument model, not a derived product.
#: ``B_S50`` carries the x1.15 red-tail correction -- see that repo's
#: ``paper/scripts/build_bayer_passbands.py``.
PASSBANDS = Path(__file__).resolve().parent / "bayer_passbands.ecsv"

#: The table ``--validate`` checks against by default: the field the S50 system was
#: derived on, from the sibling repo. Pass ``--validate <path>`` for anything else.
DEFAULT_REFERENCE = (Path(__file__).resolve().parents[2] / "photometry-adventures-2"
                     / "paper" / "data" / "gaia_20Cep_bayer.ecsv")

#: Band order is :data:`seestar_photometry.frames.BANDS`, and the column names are per
#: unit model: the S30, S30 Pro and S50 Pro use different sensors and get their own.
BANDS = ("R", "G", "B")
S50_COLUMNS = tuple(f"{band.lower()}_s50_mag{suffix}"
                    for band in BANDS for suffix in ("", "_err"))

#: Characters that turn a row of the bulk ECSV into plain wide CSV. ``flux`` and
#: ``flux_error`` arrive as quoted bracketed lists; dropping the syntax leaves
#: ``4 + 2 * 343`` numeric fields that pyarrow's C++ reader parses natively. Measured on
#: a 30 MB file: 1.4 s against 5.7 s for ``Table.read``, and bit-identical output.
_LIST_SYNTAX = b'"[]'


# --- the bulk listing -------------------------------------------------------------------

def enumerate_files(work=None, refresh=False):
    """Every bulk file, as ``(name, hp8_lo, hp8_hi, bytes)`` sorted by ``hp8_lo``.

    Cached in ``work`` after the first call. Fetching it is two seconds and four
    requests, but the two stages must agree on the same list, and a listing that changed
    between them would silently leave a hole in the sky.
    """
    cache = Path(work) / LISTING if work is not None else None
    if cache is not None and cache.exists() and not refresh:
        return [tuple(row) for row in json.loads(cache.read_text(encoding="utf-8"))]

    entries, marker = [], ""
    while True:
        url = f"{BUCKET}?prefix={LISTING_PREFIX}" + (f"&marker={marker}" if marker else "")
        with urllib.request.urlopen(url, timeout=120) as response:
            xml = response.read().decode()
        keys = re.findall(r"<Key>([^<]+\.csv\.gz)</Key>", xml)
        sizes = [int(s) for s in re.findall(r"<Size>(\d+)</Size>", xml)]
        if not keys:
            break
        for key, size in zip(keys, sizes):
            name = key.rsplit("/", 1)[1]
            lo, hi = (int(x) for x in name.split("_")[1].split(".")[0].split("-"))
            entries.append((name, lo, hi, size))
        if "<IsTruncated>true</IsTruncated>" not in xml:
            break
        marker = keys[-1]

    entries.sort(key=lambda row: row[1])
    # The product is a contiguous tiling of the sky and the join relies on it -- a gap
    # would show up as a patch of masked magnitudes rather than as an error.
    for (_, _, hi, _), (_, lo, _, _) in zip(entries, entries[1:]):
        if lo != hi + 1:
            raise RuntimeError(f"gap in the bulk listing: {hi} -> {lo}")
    if entries[0][1] != 0 or entries[-1][2] != 12 * 4 ** HPX8_LEVEL - 1:
        raise RuntimeError(f"listing does not cover the sky: {entries[0][1]}-{entries[-1][2]}")

    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(entries), encoding="utf-8")
    return entries


def files_for_parts(entries, lo_part, hi_part):
    """The listing rows overlapping HEALPix parts ``lo_part..hi_part`` inclusive."""
    lo = lo_part * HP8_PER_PART
    hi = (hi_part + 1) * HP8_PER_PART - 1
    return [row for row in entries if row[1] <= hi and row[2] >= lo]


def local_path(name, cache_dir):
    return Path(cache_dir) / name


def fetch_one(name, size, cache_dir, quiet=False, attempts=FETCH_ATTEMPTS):
    """Download one bulk file unless it is already there at the right size.

    Size rather than a checksum: the product ships no per-file digest, and a truncated
    transfer is the only realistic corruption over HTTPS. It is also the failure that
    actually happens -- one file arrived 12.5 MB short of its 31.8 MB, 28 GB into the
    first full fetch, and then downloaded cleanly three times in a row a minute later. So
    a short read is retried rather than raised, on the same reasoning as
    :func:`pipeline.solve_all`: retries belong on network stages and nowhere else.
    """
    path = local_path(name, cache_dir)
    if path.exists() and path.stat().st_size == size:
        return path, 0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    for attempt in range(1, attempts + 1):
        t0 = time.time()
        try:
            with urllib.request.urlopen(CDN + name, timeout=600) as response, \
                    open(tmp, "wb") as fh:
                while block := response.read(1 << 22):
                    fh.write(block)
            got = tmp.stat().st_size
            if got != size:
                raise OSError(f"got {got} bytes, listing says {size}")
        except (OSError, TimeoutError) as exc:
            # Half a file left behind would be skipped by the next run as "already
            # there" if it happened to be the right length, and is useless otherwise.
            tmp.unlink(missing_ok=True)
            if attempt == attempts:
                raise RuntimeError(
                    f"{name}: {type(exc).__name__}: {exc} (after {attempts} attempts)"
                ) from exc
            time.sleep(2 ** attempt)
            continue
        tmp.replace(path)
        if not quiet:
            print(f"  {name}: {size/1e6:5.0f} MB in {time.time()-t0:5.1f} s "
                  f"({size/1e6/max(time.time()-t0, 1e-3):4.1f} MB/s)", flush=True)
        return path, size


def fetch_all(entries, cache_dir, workers=16, quiet=False):
    """Download every listed file. Idempotent, so an interrupted fetch resumes.

    One file that fails all its attempts does not sink the batch -- the failures are
    collected and raised at the end, by which time the other three thousand are on disk
    and a re-run has only the stragglers left to do.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    todo = [row for row in entries
            if not (local_path(row[0], cache_dir).exists()
                    and local_path(row[0], cache_dir).stat().st_size == row[3])]
    total = sum(row[3] for row in todo)
    print(f"  {len(todo)} of {len(entries)} files to fetch ({total/1e9:.1f} GB)",
          flush=True)
    t0, done, failures = time.time(), 0, []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(fetch_one, name, size, cache_dir, quiet)
                   for name, _, _, size in todo]
        for n, future in enumerate(as_completed(futures), 1):
            try:
                done += future.result()[1]
            except RuntimeError as exc:
                failures.append(str(exc))
                print(f"  FAILED {exc}", flush=True)
            if n % 50 == 0:
                rate = done / 1e6 / max(time.time() - t0, 1e-3)
                left = (total - done) / 1e6 / max(rate, 1e-3) / 60
                print(f"  {n}/{len(todo)} files, {done/1e9:.1f} GB, "
                      f"{rate:.1f} MB/s, ~{left:.0f} min left", flush=True)
    if failures:
        raise RuntimeError(f"{len(failures)} of {len(todo)} files could not be fetched; "
                           "re-run to retry only those")
    return total


# --- parsing and integration -------------------------------------------------------------

def body_offset(raw):
    """Byte offset of the first data row in a decompressed bulk file.

    Skips Gaia's ``#`` provenance block and the column-name line after it. Counted per
    file rather than hardcoded: the block's length is not part of any contract, and
    pyarrow's CSV reader has no notion of a comment character.
    """
    newline = b"\n"
    start = 0
    while raw.startswith(b"#", start):
        start = raw.index(newline, start) + 1
    if not start:
        raise ValueError("no comment block; this is not a Gaia bulk file")
    return raw.index(newline, start) + 1        # and past the column-name line


def parse_file(path, wanted=None):
    """``(source_id, flux, flux_error)`` from one bulk file, as float32 arrays.

    ``wanted`` is a **sorted** array of ``source_id`` to keep; ``None`` keeps everything.
    Filtering happens after the parse because the rows are only identifiable once
    parsed, but before anything is copied into the output arrays.

    Bit-identical to ``Table.read(path, format="ascii.ecsv")`` and four times faster --
    see :data:`_LIST_SYNTAX` for why.
    """
    import pyarrow as pa
    import pyarrow.csv as pacsv

    with gzip.open(path, "rb") as fh:
        raw = fh.read()

    body = raw[body_offset(raw):].translate(bytes.maketrans(b"", b""), _LIST_SYNTAX)
    del raw

    n_cols = N_LEAD + 2 * N_WL
    table = pacsv.read_csv(
        io.BytesIO(body),
        read_options=pacsv.ReadOptions(autogenerate_column_names=True,
                                       block_size=1 << 26),
        convert_options=pacsv.ConvertOptions(
            column_types={"f0": pa.int64(),
                          **{f"f{i}": pa.float32() for i in range(N_LEAD, n_cols)}}),
    )
    if table.num_columns != n_cols:
        raise RuntimeError(f"{Path(path).name}: {table.num_columns} fields, "
                           f"expected {n_cols}")

    ids = table.column(0).to_numpy()
    keep = slice(None)
    if wanted is not None:
        if not len(wanted):
            return ids[:0], np.empty((0, N_WL), np.float32), np.empty((0, N_WL), np.float32)
        where = np.clip(np.searchsorted(wanted, ids), 0, len(wanted) - 1)
        keep = np.flatnonzero(wanted[where] == ids)
        ids = ids[keep]

    flux = np.empty((len(ids), N_WL), np.float32)
    error = np.empty((len(ids), N_WL), np.float32)
    for i in range(N_WL):
        flux[:, i] = table.column(N_LEAD + i).to_numpy()[keep]
        error[:, i] = table.column(N_LEAD + N_WL + i).to_numpy()[keep]
    return ids, flux, error


def load_passbands(path=PASSBANDS):
    """The three responses interpolated onto :data:`XP_GRID`, zero outside their range."""
    from astropy.table import Table

    bands = Table.read(str(path))
    return {band: np.interp(XP_GRID, bands["wavelength_nm"], bands[f"{band}_S50"],
                            left=0.0, right=0.0)
            for band in BANDS}


def synth_mags(flux, flux_err, wl, resp):
    """Photon-weighted AB magnitude and error through one passband.

    ``flux``, ``flux_err``: ``(n_src, n_wl)`` calibrated XP spectra in W m^-2 nm^-1.

    The zero-point convention is irrelevant downstream -- the per-frame ZP fit absorbs
    any constant -- only the passband shape matters. The error propagates the per-sample
    flux errors and ignores their correlations, so it is a mild underestimate; good
    enough for weighting a calibrator.

    Verbatim from ``photometry-adventures-2/paper/scripts/gaia_bayer_synth.py``, which is
    what the validated 20 Cep numbers were produced with.
    """
    f_ab = F_AB * C_NM / wl**2                     # AB source as f_lambda
    weight = resp * wl                             # photon counting
    signal = np.trapezoid(flux * weight, wl, axis=1)
    ref = np.trapezoid(f_ab * weight, wl)
    # A spectrum can integrate to zero or below in a band the source has no light in.
    # That is a missing magnitude, not an error: it comes out nan and is written masked.
    with np.errstate(invalid="ignore", divide="ignore"):
        mag = -2.5 * np.log10(signal / ref)
        sig = np.sqrt(np.sum((flux_err * weight * np.gradient(wl))**2, axis=1))
        mag_err = 2.5 / np.log(10) * sig / signal
    return mag, mag_err


def mags_for(flux, flux_err, responses):
    """``{column: array}`` of the six S50 columns for one block of spectra."""
    out = {}
    for band in BANDS:
        mag, mag_err = synth_mags(flux, flux_err, XP_GRID, responses[band])
        out[f"{band.lower()}_s50_mag"] = np.round(mag, 4)
        out[f"{band.lower()}_s50_mag_err"] = np.round(mag_err, 4)
    return out


# --- the synthesis stage ------------------------------------------------------------------

def part_groups(lo_part, hi_part, size):
    """``(lo, hi)`` inclusive part ranges covering ``lo_part..hi_part``."""
    return [(lo, min(lo + size - 1, hi_part))
            for lo in range(lo_part, hi_part + 1, size)]


def synth_group(group, src, out, cache_dir, entries, quiet=False):
    """Write the v2 parts for one contiguous run of HEALPix parts.

    Returns ``(parts_written, rows, rows_with_spectra)``. Parts already present in
    ``out`` are skipped, which is what makes the stage resumable at part granularity.
    """
    from astropy.table import MaskedColumn

    lo_part, hi_part = group
    pixels = [p for p in range(lo_part, hi_part + 1)
              if gaiadb.part_path(p, src).exists()
              and not gaiadb.part_path(p, out).exists()]
    if not pixels:
        return 0, 0, 0

    tables = {p: gaiadb.read_part(gaiadb.part_path(p, src)) for p in pixels}
    # One sorted id array for the whole group: parts are written sorted by source_id and
    # the identifier's high bits are the pixel, so concatenating them in pixel order is
    # already globally sorted and searchsorted needs no further work.
    wanted = np.concatenate([np.asarray(tables[p]["source_id"], np.int64) for p in pixels])
    responses = load_passbands()

    found = {}
    for name, _, _, _ in files_for_parts(entries, lo_part, hi_part):
        path = local_path(name, cache_dir)
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is missing; run --stage fetch before --stage synth")
        ids, flux, flux_err = parse_file(path, wanted)
        if not len(ids):
            continue
        block = mags_for(flux, flux_err, responses)
        pixel_of = gaiadb.hpx_of(ids)
        for pixel in np.unique(pixel_of):
            take = np.flatnonzero(pixel_of == pixel)
            entry = found.setdefault(int(pixel), {"source_id": []})
            entry["source_id"].append(ids[take])
            for column, values in block.items():
                entry.setdefault(column, []).append(values[take])

    written = rows = matched = 0
    for pixel in pixels:
        table = tables[pixel]
        ids = np.asarray(table["source_id"], np.int64)
        got = found.get(pixel)
        if got is None:
            values = {c: np.full(len(ids), np.nan) for c in S50_COLUMNS}
        else:
            have = np.concatenate(got["source_id"])
            order = np.argsort(have, kind="stable")
            have = have[order]
            where = np.clip(np.searchsorted(have, ids), 0, len(have) - 1)
            hit = have[where] == ids
            values = {}
            for column in S50_COLUMNS:
                column_values = np.concatenate(got[column])[order]
                values[column] = np.where(hit, column_values[where], np.nan)
            matched += int(hit.sum())
        for column in S50_COLUMNS:
            finite = np.isfinite(values[column])
            table[column] = MaskedColumn(np.nan_to_num(values[column]).astype("f4"),
                                         mask=~finite, name=column, unit="mag")
        gaiadb.write_part(table, out, pixel=pixel)
        written += 1
        rows += len(table)

    if not quiet:
        print(f"  parts {lo_part:5d}-{hi_part:5d}: {written:3d} written, "
              f"{rows/1e3:8.1f}k rows, {100*matched/max(rows, 1):5.1f}% with spectra",
              flush=True)
    return written, rows, matched


def _worker(args):
    group, src, out, cache_dir, entries = args
    return synth_group(group, src, out, cache_dir, entries)


def synth_all(src, out, cache_dir, entries, lo_part, hi_part, group_size, workers):
    """Run the synthesis stage over a range of parts.

    Reckon on ~1.2 GB of resident memory per worker: a group holds one parsed bulk file
    at a time and the largest is 122 MB gzipped, which is ~340 MB of text and two
    ``(37000, 343)`` float32 arrays on top of it.
    """
    from contextlib import nullcontext

    groups = part_groups(lo_part, hi_part, group_size)
    jobs = [(g, src, out, cache_dir, entries) for g in groups]
    print(f"  {len(groups)} groups of up to {group_size} parts, {workers} workers",
          flush=True)

    written = rows = matched = 0
    if workers == 1:
        # Threads would do here too, but a single worker is the debugging path and an
        # in-process loop keeps a traceback in one piece.
        pool, results = nullcontext(), map(_worker, jobs)
    else:
        from concurrent.futures import ProcessPoolExecutor
        pool = ProcessPoolExecutor(max_workers=workers)
        results = pool.map(_worker, jobs)
    with pool:
        for w, r, m in results:
            written += w
            rows += r
            matched += m
    return written, rows, matched


# --- validation -----------------------------------------------------------------------------

def validate(work, reference_path, quiet=False):
    """Reproduce the 20 Cep synthetic magnitudes and compare against that repo's file.

    Identical inputs and an identical formula, so the two should agree to the rounding:
    anything above 0.0001 mag means the passbands or the integration have drifted.
    """
    from astropy.table import Table

    cache_dir = Path(work) / "files"
    reference = Table.read(str(reference_path))
    ids = np.sort(np.asarray(reference["source_id"], np.int64))

    entries = enumerate_files(work)
    responses = load_passbands()
    hp8 = np.unique(ids >> HPX8_SHIFT)
    needed = [row for row in entries
              if np.any((hp8 >= row[1]) & (hp8 <= row[2]))]
    print(f"  {len(needed)} bulk files cover the field "
          f"({sum(r[3] for r in needed)/1e6:.0f} MB)", flush=True)

    got_ids, blocks = [], []
    for name, _, _, size in needed:
        path, _ = fetch_one(name, size, cache_dir, quiet)
        found, flux, flux_err = parse_file(path, ids)
        if not len(found):
            continue
        got_ids.append(found)
        blocks.append(mags_for(flux, flux_err, responses))
    got_ids = np.concatenate(got_ids)
    mine = Table({"source_id": got_ids,
                  **{c: np.concatenate([b[c] for b in blocks]) for c in S50_COLUMNS}})

    merged = Table(reference["source_id", "bayer_r_mag", "bayer_r_mag_err",
                             "bayer_g_mag", "bayer_g_mag_err",
                             "bayer_b_mag", "bayer_b_mag_err"])
    from astropy.table import join
    merged = join(merged, mine, keys="source_id", join_type="inner")
    print(f"  {len(merged)} of {len(reference)} reference rows matched", flush=True)

    worst = 0.0
    for band in BANDS:
        low = band.lower()
        for suffix in ("", "_err"):
            a = np.asarray(merged[f"bayer_{low}_mag{suffix}"], float)
            b = np.asarray(merged[f"{low}_s50_mag{suffix}"], float)
            diff = np.abs(a - b)
            worst = max(worst, float(np.nanmax(diff)))
            print(f"    {low}_s50_mag{suffix:<4}: max |d| {np.nanmax(diff):.6f}  "
                  f"median |d| {np.nanmedian(diff):.6f}", flush=True)
    ok = worst <= 1e-4
    print(f"  {'PASS' if ok else 'FAIL'}: worst disagreement {worst:.6f} mag "
          "(tolerance 0.0001)", flush=True)
    return ok


# --- entry point --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", type=Path, required=True,
                    help="scratch root; bulk files go in <work>/files")
    ap.add_argument("--src", type=Path, help="the v1 dataset to add columns to")
    ap.add_argument("--out", type=Path, help="the v2 dataset to write")
    ap.add_argument("--parts", default=f"0-{gaiadb.N_PARTS - 1}",
                    help="HEALPix part range LO-HI, inclusive")
    ap.add_argument("--group", type=int, default=32,
                    help="parts per synthesis task; larger re-parses fewer files")
    ap.add_argument("--workers", type=int, default=max((os.cpu_count() or 2) // 2, 1),
                    help="synthesis processes; reckon on ~1.2 GB resident each")
    # ESA throttles a single connection to ~4 MB/s, so the fetch is entirely about
    # concurrency. Measured over a wired link: 4.4 MB/s on one stream, 50 on 16, and
    # flat above that (46 at 24, 50 at 32, 44 at 64). 16 is where it plateaus.
    ap.add_argument("--download-workers", type=int, default=16)
    ap.add_argument("--stage", choices=("list", "fetch", "synth", "all"), default="all")
    ap.add_argument("--validate", nargs="?", const=str(DEFAULT_REFERENCE),
                    help="reproduce the magnitudes in this reference table and stop")
    ap.add_argument("--refresh-listing", action="store_true")
    args = ap.parse_args()

    cache_dir = args.work / "files"
    args.work.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    entries = enumerate_files(args.work, refresh=args.refresh_listing)
    total = sum(row[3] for row in entries)
    print(f"listing: {len(entries)} files, {total/1e9:.1f} GB, "
          f"hp8 {entries[0][1]}-{entries[-1][2]}", flush=True)

    if args.validate:
        sys.exit(0 if validate(args.work, args.validate) else 1)
    if args.stage == "list":
        have = sum(p.stat().st_size for p in cache_dir.glob("*.csv.gz")) if \
            cache_dir.exists() else 0
        print(f"cached locally: {have/1e9:.1f} GB in {args.work}")
        return

    lo_part, hi_part = (int(x) for x in args.parts.split("-"))
    if args.stage in ("fetch", "all"):
        wanted = files_for_parts(entries, lo_part, hi_part)
        print(f"stage fetch: {len(wanted)} files "
              f"({sum(r[3] for r in wanted)/1e9:.1f} GB) -> {cache_dir}", flush=True)
        fetch_all(wanted, cache_dir, workers=args.download_workers)
        print(f"stage fetch done in {(time.time()-t_start)/60:.0f} min\n", flush=True)

    if args.stage in ("synth", "all"):
        if args.src is None or args.out is None:
            ap.error("--src and --out are required for the synth stage")
        print(f"stage synth: parts {lo_part}-{hi_part} of {args.src} -> {args.out}",
              flush=True)
        t0 = time.time()
        written, rows, matched = synth_all(args.src, args.out, cache_dir, entries,
                                           lo_part, hi_part, args.group, args.workers)
        gaiadb.finalise_manifest(args.out)
        print(f"\nstage synth done: {written} parts, {rows/1e6:.1f}M rows, "
              f"{100*matched/max(rows, 1):.1f}% with spectra, "
              f"in {(time.time()-t0)/60:.0f} min", flush=True)
        size = sum(p.stat().st_size for p in args.out.rglob("*.parquet"))
        per_row = f", {size/rows:.1f} B/row" if rows else ""
        print(f"dataset: {size/1e9:.2f} GB{per_row}")
    print(f"total: {(time.time()-t_start)/60:.0f} min")


if __name__ == "__main__":
    main()
