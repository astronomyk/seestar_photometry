"""The offline Gaia reference catalogue: a local copy of what the TAP query returns.

:mod:`catalogs` builds its per-field catalogue from a Gaia TAP query. That query is the
least reliable step in the whole package -- it truncates under load, which is why
:func:`catalogs.fetch_gaia_mosaic` retries four times and refuses to cache a partial
result. This module removes it: given a local copy of the catalogue, building a field is
a file read, and the measurement path touches the network exactly never.

It also makes :func:`astrometry.solve_local` possible. A plate solve against a local
catalogue needs the catalogue to be *there* before the frame is solved, which a
per-field TAP query cannot provide -- it is the same chicken-and-egg the on-board WCS
was supposed to break.

**This is an opt-in download.** Nothing in a default install touches it, and every
function here fails with an actionable message rather than silently fetching gigabytes.

What is in it
-------------

Rows are the Gaia DR3 sources carrying synthetic photometry
(``gaiadr3.synthetic_photometry_gspc``, itself limited to G < 17.65), cut at
**V < 17.5**. Restricting to sources that have synthetic V means every row is a usable
calibrator: there is no dead weight, and nothing is masked in the columns that matter.
Cutting on V rather than G costs only the ~10-15% of rows that are red enough for V to
fall below G, and buys a faint limit that means something in the science band.

Columns are the ones the package reads today plus the ones it is about to want, and no
more -- every extra float32 is half a gigabyte over the full sky. See ``COLUMNS``.

Layout
------

A Hive-partitioned Parquet dataset, one file per HEALPix level-5 pixel (12288 of them,
under a megabyte each), rows sorted by ``source_id`` within a part::

    gaia-seestar-dr3-v1/
      MANIFEST.json
      hpx5=00000/part.parquet
      ...
      hpx5=12287/part.parquet

Partitioning at all is what makes a *partial* copy worth having: a user observing three
targets fetches tens of MB rather than the whole sky. Level 5 rather than something
coarser because of how the cone query has to work -- see :func:`parts_in_cone`, where
the pixel size sets a floor on how much sky a query must read. Measured over-read for a
1.5 degree cone: 12.4x the ideal at level 4, 5.1x at level 5, 2.7x at level 6. Level 6
would be better still but quadruples the file count to 49152 and shrinks the parts to
~150 kB, which is below the size Parquet is efficient at.

Deriving the partition needs no HEALPix arithmetic at all. A Gaia ``source_id`` encodes
its level-12 nested HEALPix index in its high bits, so ``source_id >> 49`` *is* the
level-5 pixel (see :func:`hpx_of`). Only the query side needs a HEALPix library, and
only to look up where the pixels are.

Where it lives, in order of precedence:

``SEESTAR_GAIA_DATA``
    An explicit directory. Point this at a pre-populated path on an offline machine.
``XDG_CACHE_HOME`` / ``LOCALAPPDATA``
    The platform cache location, if set.
otherwise
    ``~/.cache/seestar-photometry``.

The sibling of :func:`examples.data_dir`, with its own override so the 18 MB of example
data and a multi-GB catalogue can live on different disks -- which at this size people
do want.
"""

import json
import os
from pathlib import Path

import numpy as np

#: Gaia data release the dataset is built from. Part of the directory name because DR4
#: changes every value here -- a bare version number would let two releases collide.
DATA_RELEASE = "dr3"

#: Bumped when the schema or the row selection changes, not when Gaia does.
DATA_VERSION = "v1"

#: Directory (and remote path component) the dataset lives under.
DATASET = f"gaia-seestar-{DATA_RELEASE}-{DATA_VERSION}"

#: Where to fetch parts from. Overridden by ``SEESTAR_GAIA_URL``, which is how you point
#: at a staging host or a ``file://`` copy without editing code. Shares a host with
#: :data:`astap.BASE_URL`; see ``tools/HOSTING.md``.
BASE_URL = f"https://crowdsky.univie.ac.at/seestar_assets/{DATASET}"

#: HEALPix order of the partitioning. See the module docstring for why 5.
HPX_LEVEL = 5

#: Number of parts a complete dataset has: ``12 * 4 ** HPX_LEVEL``.
N_PARTS = 12 * 4 ** HPX_LEVEL

#: Bits to shift a ``source_id`` right to get its HEALPix pixel. Gaia packs the level-12
#: nested index above 35 bits of running number, and level 12 -> level 5 is a further
#: ``2 * (12 - 5)`` bits.
_HPX_SHIFT = 35 + 2 * (12 - HPX_LEVEL)

#: Greatest angular distance from a pixel centre to any point inside that pixel, in
#: degrees, at :data:`HPX_LEVEL`. Measured over all 12288 pixel boundaries rather than
#: derived, because HEALPix pixels are not circles and the polar ones are the worst
#: case: 1.903 deg against a 1.833 deg nominal half-resolution.
#:
#: This is what makes :func:`parts_in_cone` provably complete, so it is rounded *up*.
HPX_MAX_RADIUS = 1.91

#: Pixel centres, computed once on first use. 12288 of them, so not worth recomputing.
_CENTRES = None

#: Epoch of the catalogue positions. Gaia DR3 is J2016.0 for every source.
REF_EPOCH = 2016.0

#: Faint limit in synthetic V. See the module docstring.
V_LIMIT = 17.5

#: Filename of the per-part checksum manifest.
MANIFEST = "MANIFEST.json"

#: Columns carried, in on-disk order. ``cone`` returns exactly these.
#:
#: The first nine are what the package reads today and match what
#: :func:`catalogs.fetch_gaia_mosaic` produces. The rest are cheap additions whose
#: consumers are still to be written: ``c_star`` and ``ipd_frac_multi_peak`` flag blends
#: (at 2.4 arcsec/px with a ~10 arcsec FWHM, blended comparisons are the norm, and this
#: is Gaia telling you which ones), ``ruwe`` and ``non_single_star`` flag unresolved
#: binaries, and ``teff_gspphot`` gives a stellar radius for the transit use case.
#:
#: ``v_jkc_flag`` is Gaia's own statement that a source's G and BP-RP fall inside the
#: range over which the synthetic V was validated. It guards the one column everything
#: else depends on, and as a boolean it compresses to nothing, so it is carried even
#: though nothing reads it yet.
COLUMNS = (
    "source_id", "ra", "dec",
    "phot_g_mean_mag", "bp_rp", "phot_variable_flag",
    "v_jkc_mag", "b_jkc_mag", "r_jkc_mag", "v_jkc_flag",
    "pmra", "pmdec",
    "c_star", "ruwe", "teff_gspphot",
    "ipd_frac_multi_peak", "non_single_star",
    "duplicated_source", "in_galaxy_candidates", "has_epoch_photometry",
)

#: Storage dtype per column. ``"str"`` is dictionary-encoded on disk.
DTYPES = {
    "source_id": "int64",
    "ra": "float64",
    "dec": "float64",
    "phot_g_mean_mag": "float32",
    "bp_rp": "float32",
    "phot_variable_flag": "str",
    "v_jkc_mag": "float32",
    "b_jkc_mag": "float32",
    "r_jkc_mag": "float32",
    "v_jkc_flag": "bool",
    "pmra": "float32",
    "pmdec": "float32",
    "c_star": "float32",
    "ruwe": "float32",
    "teff_gspphot": "float32",
    "ipd_frac_multi_peak": "uint8",
    "non_single_star": "uint8",
    "duplicated_source": "bool",
    "in_galaxy_candidates": "bool",
    "has_epoch_photometry": "bool",
}

#: Units, so a table read from here is indistinguishable from a TAP result.
UNITS = {
    "ra": "deg", "dec": "deg",
    "pmra": "mas / yr", "pmdec": "mas / yr",
    "phot_g_mean_mag": "mag", "bp_rp": "mag",
    "v_jkc_mag": "mag", "b_jkc_mag": "mag", "r_jkc_mag": "mag",
    "teff_gspphot": "K",
}

#: Columns that are never null in the source data, so never masked on the way out.
_ALWAYS_PRESENT = ("source_id", "ra", "dec")

#: Value written under a null before the mask is applied. Only ever seen through the
#: mask, but it must be castable to the column's dtype -- NaN is not an integer and
#: casts to ``True`` as a bool, which is exactly the kind of silent wrongness the mask
#: exists to prevent.
_FILL = {"str": "", "bool": False}

#: numpy dtype to build an all-missing column with, where the schema's own is awkward
#: to allocate empty. Only used for columns a dataset does not carry.
_EMPTY_DTYPE = {"str": "U16"}


def hpx_of(source_id):
    """Level-:data:`HPX_LEVEL` HEALPix pixel of one or many Gaia ``source_id``.

    Pure integer arithmetic: Gaia packs the nested level-12 index into the high bits of
    the identifier, so no HEALPix library and no trigonometry is involved. This is what
    makes the dataset cheap to partition at build time.
    """
    return np.asarray(source_id, dtype=np.int64) >> _HPX_SHIFT


def separation_deg(ra0, dec0, ra, dec):
    """Great-circle separation in degrees between one point and an array of points.

    Haversine rather than the dot-product form, which loses precision at exactly the
    small separations a cone cut cares about.
    """
    lat0, lat = np.radians(dec0), np.radians(np.asarray(dec, dtype=float))
    dlat = (lat - lat0) / 2.0
    dlon = (np.radians(np.asarray(ra, dtype=float)) - np.radians(ra0)) / 2.0
    a = np.sin(dlat) ** 2 + np.cos(lat0) * np.cos(lat) * np.sin(dlon) ** 2
    return np.degrees(2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))))


# --- where the data lives --------------------------------------------------------------

def cache_dir():
    """Root directory the catalogue is cached under. Downloads nothing.

    Honours ``SEESTAR_GAIA_DATA`` first, then the platform cache location.
    """
    override = os.environ.get("SEESTAR_GAIA_DATA")
    if override:
        return Path(override).expanduser()
    base = os.environ.get("XDG_CACHE_HOME") or os.environ.get("LOCALAPPDATA")
    root = Path(base).expanduser() if base else Path.home() / ".cache"
    return root / "seestar-photometry"


def dataset_dir():
    """Directory of this release's dataset. Downloads nothing."""
    return cache_dir() / DATASET


def part_path(pixel, directory=None):
    """Path of one HEALPix part, present or not."""
    directory = Path(directory) if directory is not None else dataset_dir()
    return directory / f"hpx{HPX_LEVEL}={int(pixel):04d}" / "part.parquet"


def _pixel_centres():
    """``(lon, lat)`` of every HEALPix pixel centre, in degrees. Computed once."""
    global _CENTRES
    if _CENTRES is None:
        import astropy.units as u
        from astropy_healpix import HEALPix

        hp = HEALPix(nside=2 ** HPX_LEVEL, order="nested", frame="icrs")
        lon, lat = hp.healpix_to_lonlat(np.arange(hp.npix))
        _CENTRES = (np.asarray(lon.to_value(u.deg), dtype=float),
                    np.asarray(lat.to_value(u.deg), dtype=float))
    return _CENTRES


def parts_in_cone(center, radius_deg):
    """HEALPix pixels overlapping a cone, as a sorted integer array.

    Every pixel whose *centre* lies within ``radius_deg + HPX_MAX_RADIUS`` of the cone
    centre. That is complete by construction: if any point of a pixel is inside the
    cone, that point is within ``radius_deg`` of the cone centre and within
    :data:`HPX_MAX_RADIUS` of its own pixel centre, so the pixel centre is inside the
    padded radius and the pixel is returned. It also covers a cone smaller than a pixel,
    since a pixel centre is never further than the padding from a point it contains.

    ``astropy_healpix.cone_search_lonlat`` is deliberately not used here: it returns
    fewer pixels than this, and under-selection silently drops real sources rather than
    raising -- which is precisely the failure the frame table would report as a
    mysteriously low ``n_cal``.

    The price of completeness is reading more sky than asked for -- about 5x for a 1.5
    degree cone. The exact cut in :func:`cone` removes the surplus, so it costs I/O and
    never accuracy.
    """
    lon, lat = _pixel_centres()
    reach = float(radius_deg) + HPX_MAX_RADIUS
    return np.flatnonzero(separation_deg(center[0], center[1], lon, lat) <= reach)


def covers(center, radius_deg, directory=None):
    """Whether a cone can be answered completely from disk. Downloads nothing.

    Every part the cone reaches must be present -- a missing one might hold sources
    inside the cone, and there is no way to know without reading it.

    A part the *dataset* does not carry is not counted as missing. That distinction
    matters both ways round. A regionally-built dataset has a manifest listing only what
    it was built from, and is complete for its own region. A full dataset that has been
    partially fetched has the whole manifest, so a part listed there but absent locally
    is genuinely missing and this returns False -- which is what makes
    ``catalogue_backend="auto"`` fall back to TAP for a field you have not fetched.
    """
    directory = Path(directory) if directory is not None else dataset_dir()
    meta = manifest(directory)
    published = set(meta["parts"]) if meta else None

    found = 0
    for pixel in parts_in_cone(center, radius_deg):
        if published is not None and str(int(pixel)) not in published:
            continue
        if not part_path(pixel, directory).exists():
            return False
        found += 1
    return found > 0


def manifest(directory=None):
    """The parsed ``MANIFEST.json``, or ``None`` if the dataset is not present."""
    directory = Path(directory) if directory is not None else dataset_dir()
    path = directory / MANIFEST
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


# --- reading ----------------------------------------------------------------------------

def _column_to_masked(chunked, name):
    """One Arrow column as an astropy column, nulls preserved as a mask.

    Load-bearing: every "did this source match" test in the package is
    ``np.ma.getmaskarray``, never ``np.isfinite``. A plain column carrying NaN reads as
    a *successful* match with a NaN magnitude, which propagates into the zero-point fit
    as a silently wrong answer rather than an error.
    """
    import pyarrow as pa
    from astropy.table import Column, MaskedColumn

    dtype = DTYPES[name]
    if pa.types.is_dictionary(chunked.type):
        chunked = chunked.cast(pa.string())
    mask = np.asarray(chunked.is_null(), dtype=bool)
    filled = chunked.fill_null(_FILL.get(dtype, 0))
    values = np.asarray(filled)
    values = values.astype(str) if dtype == "str" else values.astype(dtype)

    unit = UNITS.get(name)
    if name in _ALWAYS_PRESENT:
        return Column(values, name=name, unit=unit)
    return MaskedColumn(values, mask=mask, name=name, unit=unit)


def _to_table(arrow):
    """An Arrow table as an astropy Table with this module's schema."""
    from astropy.table import Table

    out = Table()
    for name in COLUMNS:
        out[name] = _column_to_masked(arrow.column(name), name)
    return out


def cone(center, radius_deg, directory=None, vmag_limit=None, gmag_limit=None,
         epoch=None):
    """Every catalogue source within ``radius_deg`` of ``center``, as an astropy Table.

    A drop-in for the result of :func:`catalogs.fetch_gaia_mosaic`: same column names,
    same units, same masking convention, plus the extra columns in :data:`COLUMNS`.

    Parameters
    ----------
    center : tuple of float
        ``(ra, dec)`` in degrees.
    radius_deg : float
        Cone radius. For a project this wants to cover the whole dithered area, not one
        field of view.
    directory : path-like, optional
        Dataset root. Defaults to the cache; pass one to read a dataset built locally.
    vmag_limit : float, optional
        Drop sources fainter than this in synthetic V. Sources with no V are dropped,
        since without it they cannot be calibrators.
    gmag_limit : float, optional
        Drop sources fainter than this in Gaia G. This is the cut the TAP query applies
        (``Project.gmag_limit``), so passing it makes the two backends select the same
        sources.
    epoch : float, optional
        Decimal year to propagate positions to. See :func:`propagate_epoch`.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    directory = Path(directory) if directory is not None else dataset_dir()
    pixels = parts_in_cone(center, radius_deg)
    paths = [p for p in (part_path(i, directory) for i in pixels) if p.exists()]
    if not paths:
        raise FileNotFoundError(
            f"no local Gaia parts for the cone at {tuple(center)} r={radius_deg} deg\n"
            f"  looked in {directory}\n"
            f"  needed HEALPix parts {list(pixels)}\n"
            "Fetch them with gaiadb.download(center=..., radius_deg=...), or point "
            "SEESTAR_GAIA_DATA at a directory that has them."
        )

    arrow = pa.concat_tables([pq.read_table(p, columns=list(COLUMNS)) for p in paths])
    table = _to_table(arrow)

    keep = separation_deg(center[0], center[1], table["ra"], table["dec"]) <= radius_deg
    # Filled with +inf, so a source missing the magnitude fails the cut rather than
    # passing it -- one with no synthetic V could not be calibrated against anyway.
    if vmag_limit is not None:
        keep &= np.ma.filled(table["v_jkc_mag"], np.inf) < vmag_limit
    if gmag_limit is not None:
        keep &= np.ma.filled(table["phot_g_mean_mag"], np.inf) < gmag_limit
    table = table[keep]
    table.meta["epoch"] = REF_EPOCH
    table.meta["catalogue"] = DATASET
    if epoch is not None:
        table = propagate_epoch(table, epoch)
    return table


def propagate_epoch(catalogue, epoch, ref_epoch=REF_EPOCH):
    """Move catalogue positions from Gaia's epoch to the observing epoch.

    Gaia DR3 positions are J2016.0. Ten years later a star with a 100 mas/yr proper
    motion has moved 1 arcsec -- half the default 2 arcsec match tolerance, and enough
    to lose the match entirely on the fastest movers. The TAP path cannot do this at all
    because it never fetches proper motions.

    Linear propagation on the tangent plane, which is right to well under a milliarcsec
    over a decade. ``pmra`` is Gaia's ``pmra*``, already multiplied by ``cos(dec)``, so
    it is divided back out to get a change in right ascension.
    """
    dt = float(epoch) - float(ref_epoch)
    dec = np.asarray(catalogue["dec"], dtype=float)
    pmra = np.ma.filled(catalogue["pmra"], 0.0).astype(float)
    pmdec = np.ma.filled(catalogue["pmdec"], 0.0).astype(float)

    out = catalogue.copy()
    out["dec"] = dec + pmdec * dt / 3.6e6
    # cos(dec) is floored rather than guarded: within a few mas of a pole right
    # ascension is meaningless anyway, and a division by zero here would produce a NaN
    # position that reads as a real source.
    cos_dec = np.maximum(np.cos(np.radians(dec)), 1e-8)
    out["ra"] = np.mod(np.asarray(catalogue["ra"], dtype=float)
                       + pmra * dt / 3.6e6 / cos_dec, 360.0)
    out.meta["epoch"] = float(epoch)
    return out


# --- writing ------------------------------------------------------------------------------

def _to_arrow(catalogue):
    """An astropy Table as an Arrow table in this module's schema."""
    import pyarrow as pa

    fields, arrays = [], []
    n_rows = len(catalogue)
    for name in COLUMNS:
        if name not in catalogue.colnames:
            if name in _ALWAYS_PRESENT:
                raise ValueError(f"catalogue is missing the {name!r} column")
            # A dataset may legitimately carry only some of the schema -- the columns
            # come from two Gaia tables and the expensive one is 790 GB in bulk form.
            # An absent column is written as entirely null, so it reads back masked,
            # which is exactly what "we do not have this" should look like to code that
            # tests the mask. The manifest records which are real.
            column = np.ma.masked_all(n_rows, dtype=_EMPTY_DTYPE.get(DTYPES[name],
                                                                     DTYPES[name]))
        else:
            column = catalogue[name]
        dtype = DTYPES[name]
        mask = np.ma.getmaskarray(column)
        values = np.ma.getdata(column)

        if dtype == "str":
            array = pa.array(np.asarray(values, dtype=str), type=pa.string(), mask=mask)
            array = array.dictionary_encode()
        elif dtype == "bool":
            array = pa.array(np.where(mask, False, values).astype(bool), mask=mask)
        elif np.issubdtype(np.dtype(dtype), np.integer):
            values = np.asarray(values)
            if np.issubdtype(values.dtype, np.floating):
                # A float input reaching an integer column carries nulls as NaN, which
                # has no integer representation -- numpy's cast of it is
                # platform-defined garbage rather than an error, so zero them first.
                #
                # Only for a float input. Laundering an integer column through float64
                # loses every bit below 2**53, and a Gaia `source_id` is ~1e18: the
                # identifier would come back subtly wrong and match nothing.
                values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
            array = pa.array(np.where(mask, 0, values).astype(dtype), mask=mask)
        else:
            array = pa.array(np.asarray(values, dtype=dtype), mask=mask)

        fields.append(pa.field(name, array.type))
        arrays.append(array)
    return pa.Table.from_arrays(arrays, schema=pa.schema(fields))


def sha256_of(path, chunk=1 << 20):
    """SHA-256 of a file, read in chunks.

    Chunked because these parts are large; ``examples`` hashes whole archives in memory,
    which does not scale to a dataset this size.
    """
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def write_dataset(catalogue, directory, compression="zstd"):
    """Write a catalogue Table out as a partitioned dataset, and return its manifest.

    Used by the full-sky build tool and by the tests, which is why it takes a Table
    rather than talking to Gaia itself -- the source of the rows is not this module's
    problem.

    Rows are sorted by ``source_id`` within each part. That is a spatial sort for free:
    the identifier's high bits are the level-12 HEALPix index, so sorting on it also
    clusters neighbours together and lets the float columns delta-encode well.

    Writing a *subset* of the sky is expected and produces a valid dataset -- the
    manifest records only the parts that exist, and :func:`covers` answers accordingly.
    """
    import pyarrow.parquet as pq

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    order = np.argsort(np.asarray(catalogue["source_id"], dtype=np.int64), kind="stable")
    catalogue = catalogue[order]
    pixels = hpx_of(catalogue["source_id"])
    for pixel in np.unique(pixels):
        write_part(catalogue[pixels == pixel], directory, compression=compression)
    return finalise_manifest(directory)


def write_part(catalogue, directory, pixel=None, compression="zstd"):
    """Write the rows of one HEALPix part, and return its path.

    Split out from :func:`write_dataset` so an all-sky build can go part by part and be
    killed and resumed: 12288 parts is many hours of work, and redoing it because a
    laptop slept is not acceptable.
    """
    import pyarrow.parquet as pq

    pixels = np.unique(hpx_of(catalogue["source_id"]))
    if pixel is None:
        if len(pixels) != 1:
            raise ValueError(
                f"write_part takes rows from one HEALPix pixel, got {len(pixels)}"
            )
        pixel = int(pixels[0])
    elif len(pixels) and pixels[0] != pixel:
        raise ValueError(f"rows belong to pixel {pixels[0]}, not {pixel}")

    arrow = _to_arrow(catalogue)
    path = part_path(pixel, directory)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".part")
    pq.write_table(
        arrow, tmp, compression=compression,
        use_byte_stream_split=[n for n in COLUMNS if DTYPES[n].startswith("float")],
        use_dictionary=["phot_variable_flag"],
    )
    # Renamed into place only once complete, so an interrupted build leaves no
    # half-written part that a resumed run would skip as already done.
    tmp.replace(path)
    return path


def finalise_manifest(directory):
    """Index whatever parts are on disk, and write the manifest.

    Separate from writing so a build can be resumed and then indexed once at the end.
    Which columns actually hold values is read from each part's Parquet null counts
    rather than tracked by the caller -- the statistics are in the footer, so it costs
    metadata reads rather than a pass over the data.
    """
    import pyarrow.parquet as pq

    directory = Path(directory)
    parts, rows, populated = {}, 0, set()
    for path in sorted(directory.glob(f"hpx{HPX_LEVEL}=*/part.parquet")):
        pixel = int(path.parent.name.split("=", 1)[1])
        meta = pq.ParquetFile(path).metadata
        parts[str(pixel)] = {
            "sha256": sha256_of(path),
            "bytes": path.stat().st_size,
            "rows": int(meta.num_rows),
        }
        rows += int(meta.num_rows)
        schema = meta.schema.to_arrow_schema()
        for group in range(meta.num_row_groups):
            for i, name in enumerate(schema.names):
                stats = meta.row_group(group).column(i).statistics
                if stats is not None and stats.null_count < meta.row_group(group).num_rows:
                    populated.add(name)

    result = {
        "dataset": DATASET,
        "release": DATA_RELEASE,
        "version": DATA_VERSION,
        "hpx_level": HPX_LEVEL,
        "ref_epoch": REF_EPOCH,
        "v_limit": V_LIMIT,
        "columns": list(COLUMNS),
        # Which of them hold real values. The rest are present in the schema and
        # entirely masked, so callers see a stable set of columns either way.
        "columns_present": [c for c in COLUMNS if c in populated],
        "rows": rows,
        "parts": parts,
    }
    (directory / MANIFEST).write_text(json.dumps(result, indent=1), encoding="utf-8")
    return result


# --- fetching -----------------------------------------------------------------------------

def _base_url():
    return os.environ.get("SEESTAR_GAIA_URL", BASE_URL).rstrip("/")


def _fetch(url, dest, expect_sha256=None, chunk=1 << 20):
    """Download one file to ``dest`` via a temporary name, verifying the checksum."""
    import urllib.error
    import urllib.request

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=120) as response, open(tmp, "wb") as fh:
            while block := response.read(chunk):
                fh.write(block)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"could not download {url}\n"
            f"  {type(exc).__name__}: {exc}\n"
            "If this machine is offline, copy the dataset across and point "
            "SEESTAR_GAIA_DATA at it."
        ) from exc

    if expect_sha256 is not None:
        digest = sha256_of(tmp)
        if digest != expect_sha256:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(
                f"checksum mismatch for {url}\n"
                f"  expected {expect_sha256}\n  got      {digest}"
            )
    tmp.replace(dest)
    return dest


def download(center=None, radius_deg=None, force=False, directory=None, quiet=False):
    """Fetch the catalogue, or just the parts covering one region.

    Idempotent: parts already present and checksum-listed are skipped, so an interrupted
    fetch resumes by being run again. Every part is verified against the manifest before
    it is moved into place, so a truncated download cannot leave a file that later reads
    as a mysteriously empty patch of sky.

    Pass ``center`` and ``radius_deg`` to fetch only the sky a project needs -- a few
    tens of MB instead of the whole catalogue. Omit both for everything.
    """
    directory = Path(directory) if directory is not None else dataset_dir()
    directory.mkdir(parents=True, exist_ok=True)

    # Always re-fetched, never reused from disk. The manifest is the index of what the
    # server carries, and that grows: a dataset published for one region and later
    # extended would otherwise be invisible to anyone who fetched the early manifest,
    # and they would be told their field is not carried when it now is. A few kB per
    # call is a cheap price for not having a stale index.
    local = directory / MANIFEST
    _fetch(f"{_base_url()}/{MANIFEST}", local)
    meta = json.loads(local.read_text(encoding="utf-8"))

    if center is not None and radius_deg is not None:
        wanted = [str(int(p)) for p in parts_in_cone(center, radius_deg)]
    else:
        wanted = sorted(meta["parts"], key=int)

    todo = []
    for pixel in wanted:
        entry = meta["parts"].get(pixel)
        if entry is None:
            continue  # a region the dataset does not cover, e.g. a partial build
        if force or not part_path(pixel, directory).exists():
            todo.append((pixel, entry))

    if not todo:
        return directory
    total = sum(e["bytes"] for _, e in todo) / 1e6
    if not quiet:
        print(f"[seestar-photometry] fetching {len(todo)} Gaia parts "
              f"(~{total:.0f} MB) -> {directory}", flush=True)
    for pixel, entry in todo:
        _fetch(f"{_base_url()}/hpx{HPX_LEVEL}={int(pixel):04d}/part.parquet",
               part_path(pixel, directory), expect_sha256=entry["sha256"])
    if not quiet:
        print(f"[seestar-photometry] Gaia catalogue ready in {directory}", flush=True)
    return directory
