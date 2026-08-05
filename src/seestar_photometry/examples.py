"""Real Seestar data, downloaded on demand, for the documentation and the tests.

Every example in the docs runs on these files, and so do the real-data tests. That is
deliberate: an example that cannot be executed rots, and a test on synthetic data alone
cannot catch a wrong assumption about what the instrument actually writes.

The frames are **1000x1000 cutouts** of real MW Cam observations, gzipped, with the WCS
kept valid (the cutout origin is subtracted from CRPIX). Solved ``.wcs`` sidecars and a
trimmed Gaia DR3 table come with them, so cross-matching, calibration and light curves all
work with **no plate solver installed**.

>>> from seestar_photometry import examples, photometry
>>> frame = examples.stack()          # downloads ~18 MB once, then cached
>>> ext = photometry.extract_sources(frame)

**Not bundled in the wheel.** 18 MB of FITS inside a package every user installs is rude to
PyPI and to anyone who only wants the photometry. The archive is a GitHub release asset,
fetched on first use and cached, so a CI runner pays the download once per cache lifetime and
an ordinary install pays nothing.

Where the cache goes, in order of precedence:

``SEESTAR_PHOTOMETRY_DATA``
    An explicit directory. Point this at a pre-populated path to skip the download
    entirely -- the right approach for an offline or air-gapped machine.
``XDG_CACHE_HOME`` / ``LOCALAPPDATA``
    The platform cache location, if set.
otherwise
    ``~/.cache/seestar-photometry``.

What is here:

:func:`stack`
    760 s c17 stack, rms 0.016 mag -- the general-purpose frame.
:func:`stack_deep`
    1460 s stack of the same field, two hours later.
:func:`stack_saturated`
    280 s stack whose brightest star clips the sensor.
:func:`crowdsky`
    A CrowdSky multi-extension frame, plate-solved server-side.
:func:`raw_subs`
    Five consecutive 20 s raw Bayer subs, for stacking.
:func:`gaia`
    Gaia DR3 rows over the footprint, with synthetic Johnson V.

Regenerate the archive with ``tools/build_example_data.py`` (needs the full datasets).
"""

import hashlib
import os
import tarfile
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

#: Version of the example dataset. Bumped only when the files change; the cache is keyed on
#: it, so a new version never needs a stale copy invalidated by hand.
DATA_VERSION = "v1"

#: Archive name, and the release tag it hangs off.
ARCHIVE = f"example-data-{DATA_VERSION}.tar.gz"

#: Where to fetch it. A release asset rather than a file in the repository, so cloning stays
#: fast and the URL is immutable once published.
DATA_URL = (
    "https://github.com/astronomyk/seestar_photometry/releases/download/"
    f"example-data-{DATA_VERSION}/{ARCHIVE}"
)

#: SHA-256 of the archive, verified after download. A truncated or tampered file would
#: otherwise surface much later as baffling photometry rather than as a download error.
DATA_SHA256 = "0c6e38cd5f8c68e6179c5eda4ccecdd5d77df61f89897c67fb573b750cc9c813"

#: What the archive must contain.
_EXPECTED = (
    "stack_c17_15min.fits.gz", "stack_c17_15min.wcs",
    "stack_c17_30min.fits.gz", "stack_c17_30min.wcs",
    "stack_saturated.fits.gz", "stack_saturated.wcs",
    "crowdsky_mef.fits.gz", "crowdsky_mef.wcs",
    "raw_sub_1.fits.gz", "raw_sub_2.fits.gz", "raw_sub_3.fits.gz",
    "raw_sub_4.fits.gz", "raw_sub_5.fits.gz",
    "gaia_mwcam.ecsv",
)

#: Field centre of every example frame (MW Cam), in degrees. Handy as a ``Target``.
FIELD_RA, FIELD_DEC = 186.6821, 81.474

_FRAMES = {
    "stack": "stack_c17_15min",
    "stack_deep": "stack_c17_30min",
    "stack_saturated": "stack_saturated",
    "crowdsky": "crowdsky_mef",
}


def data_dir():
    """Directory the example data lives in. Downloads nothing.

    Honours ``SEESTAR_PHOTOMETRY_DATA`` first, then the platform cache location.
    """
    override = os.environ.get("SEESTAR_PHOTOMETRY_DATA")
    if override:
        return Path(override).expanduser()
    base = os.environ.get("XDG_CACHE_HOME") or os.environ.get("LOCALAPPDATA")
    root = Path(base).expanduser() if base else Path.home() / ".cache"
    return root / "seestar-photometry" / DATA_VERSION


def is_downloaded(directory=None):
    """Whether every expected file is already present."""
    directory = Path(directory) if directory else data_dir()
    return all((directory / name).exists() for name in _EXPECTED)


def download(force=False, url=None, quiet=False):
    """Fetch and unpack the example data, returning its directory.

    Idempotent: returns straight away if the files are already there. Every accessor calls
    it, so you rarely need it directly -- but calling it once up front is the polite thing to
    do in a test session or a CI step, so the download is not attributed to whichever test
    happened to run first.

    The archive is checksum-verified before being unpacked, and unpacked through a staging
    directory so an interrupted download cannot leave a half-populated cache that
    :func:`is_downloaded` would then accept as complete.
    """
    target = data_dir()
    if not force and is_downloaded(target):
        return target

    url = url or DATA_URL
    target.mkdir(parents=True, exist_ok=True)
    if not quiet:
        print(f"[seestar-photometry] fetching example data (~18 MB) -> {target}",
              flush=True)

    with tempfile.TemporaryDirectory(prefix="sp_data_") as tmp:
        archive = Path(tmp) / ARCHIVE
        try:
            with urllib.request.urlopen(url, timeout=120) as response, \
                    open(archive, "wb") as fh:
                while chunk := response.read(1 << 20):
                    fh.write(chunk)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            raise RuntimeError(
                f"could not download the example data from {url}\n"
                f"  {type(exc).__name__}: {exc}\n"
                "If this machine is offline, fetch the archive elsewhere, unpack it, and "
                "point SEESTAR_PHOTOMETRY_DATA at that directory."
            ) from exc

        digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        if digest != DATA_SHA256:
            raise RuntimeError(
                f"example data checksum mismatch for {url}\n"
                f"  expected {DATA_SHA256}\n  got      {digest}\n"
                "The download is corrupt, or the release asset was replaced."
            )

        staging = Path(tmp) / "unpacked"
        staging.mkdir()
        with tarfile.open(archive, "r:gz") as tf:
            for member in tf.getmembers():
                # The archive is flat by construction; refuse anything trying to escape.
                if member.name != Path(member.name).name or not member.isfile():
                    raise RuntimeError(f"unexpected archive member {member.name!r}")
                tf.extract(member, staging)
        for name in _EXPECTED:
            source = staging / name
            if not source.exists():
                raise RuntimeError(f"example data archive is missing {name!r}")
            source.replace(target / name)

    if not quiet:
        print(f"[seestar-photometry] example data ready in {target}", flush=True)
    return target


def path(name):
    """Path to an example file by stem, downloading the dataset if needed.

    Raises a listing rather than a bare ``FileNotFoundError``, because the usual cause is a
    typo and the useful answer is what *is* available.
    """
    directory = download()
    for candidate in (f"{name}.fits.gz", f"{name}.ecsv", f"{name}.wcs", name):
        hit = directory / candidate
        if hit.exists():
            return hit
    raise FileNotFoundError(
        f"no example file {name!r}. Available: {', '.join(available())}"
    )


def available():
    """Stems of every example file present in the cache. Downloads nothing."""
    directory = data_dir()
    if not directory.is_dir():
        return []
    return sorted({p.name.split(".")[0] for p in directory.iterdir() if p.is_file()})


def _load(key):
    from .frames import load_frame

    return load_frame(path(_FRAMES[key]))


def stack():
    """A clean 760 s c17 stack: rms 0.016 mag, 59 calibration stars, FWHM 3.8 px.

    The frame to reach for first. Native ``"cube"`` layout, with a solved sidecar.
    """
    return _load("stack")


def stack_deep():
    """A 1460 s stack of the same field, two hours later.

    Paired with :func:`stack` it gives a real depth-vs-exposure comparison and a
    two-epoch light curve -- the same stars, the same unit, a different integration.
    """
    return _load("stack_deep")


def stack_saturated():
    """A 280 s stack whose brightest star clips the sensor.

    Needed for the saturation-limit use case: ``calibration.saturation_mag`` returns
    ``nan`` on a frame where nothing saturates, so demonstrating it takes a frame that
    does. This one's bright limit is V ~ 8.1.
    """
    return _load("stack_saturated")


def crowdsky():
    """A CrowdSky frame: multi-extension layout, plate-solved server-side.

    The second FITS layout, with the different header dialect (``NIMAGES``, and
    ``EXPTIME`` meaning the *total* on-sky time) plus a ``FOOTPRINT`` plane and the
    server's own ``STAR-TAB`` catalogue. Its header WCS is trustworthy, so
    ``astrometry.lift`` applies.
    """
    return _load("crowdsky")


def raw_subs(n=5):
    """Up to five consecutive 20 s raw Bayer sub-exposures, as loaded frames.

    Same unit and night as :func:`stack`. Raw subs are single-plane mosaics, demosaiced on
    load; ``frame.bayer`` keeps the native samples. Feed the *paths* to
    :func:`stacking.stack_frame` -- see :func:`raw_sub_paths`.
    """
    from .frames import load_frame

    return [load_frame(p) for p in raw_sub_paths(n)]


def raw_sub_paths(n=5):
    """Paths of the raw subs, in time order -- what :mod:`stacking` wants."""
    directory = download()
    out = [directory / f"raw_sub_{i}.fits.gz" for i in range(1, n + 1)]
    missing = [p.name for p in out if not p.exists()]
    if missing:
        raise FileNotFoundError(f"example raw subs missing: {', '.join(missing)}")
    return out


def gaia():
    """Gaia DR3 rows over the example footprint, with synthetic Johnson-Kron-Cousins V.

    Trimmed from a 1.5 degree mosaic to what these cutouts need, so
    :func:`calibration.fit_zeropoint` works with no TAP query. Drop-in wherever the
    pipeline wants a reference catalogue.
    """
    from astropy.table import Table

    return Table.read(path("gaia_mwcam"))


def wcs(name="stack"):
    """The cached WCS of an example frame, as an :class:`astropy.wcs.WCS`.

    Equivalent to ``astrometry.load_wcs(frame)``; offered directly for examples that want
    a WCS without loading pixels.
    """
    from . import astrometry
    from .frames import load_frame

    frame = load_frame(path(_FRAMES[name]))
    solved = astrometry.load_wcs(frame)
    if solved is None:
        raise FileNotFoundError(f"no .wcs sidecar for {name!r}")
    return solved


def target():
    """A :class:`project.Target` for MW Cam, the field every example frame shows."""
    from .project import Target

    return Target("MW Cam", ra=FIELD_RA, dec=FIELD_DEC)


def project(work_dir, **kwargs):
    """A ready-to-run :class:`project.Project` over the example frames.

    The quickest way to exercise the batch stages end to end. The Gaia cache is seeded from
    the example table, so no query is attempted.

    >>> proj = examples.project(tmp_path)          # doctest: +SKIP
    >>> frames = pipeline.build_frame_table(proj)  # doctest: +SKIP
    """
    import shutil

    from .frames import LocalTree
    from .project import Project

    directory = download()
    proj = Project(
        target=target(),
        source=LocalTree(roots=[directory], patterns=("stack_*.fits.gz",)),
        work_dir=work_dir,
        **kwargs,
    )
    proj.ensure_dirs()
    if not proj.catalogue_path.exists():
        shutil.copy(path("gaia_mwcam"), proj.catalogue_path)
    return proj
