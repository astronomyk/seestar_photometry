"""Finding, and if necessary fetching, the ASTAP plate solver.

ASTAP is the default solver (:func:`astrometry.solve_astap`) and it is the right
default: local, offline, tried and tested, and several times faster than anything else
here. Its only drawback was installation -- a binary and a star database to download by
hand, from a path this package then had to be told about. This module removes that.

Nothing here is needed if ASTAP is already installed. :func:`executable` looks for one
in the obvious places first, so ``apt install astap-cli`` or a manual install is found
without configuration. :func:`download` is the fallback, for machines where installing
system packages is not an option::

    from seestar_photometry import astap

    astap.download()        # ~100 MB, once
    astap.executable()      # -> the path solve_astap will use

Why a download and not a dependency: there is no PyPI package that bundles ASTAP.
``astapy`` wraps the CLI but still expects you to install it by hand, which is the
problem rather than the solution. ASTAP is `MPL-2.0
<https://github.com/han-k59/astap>`_, so redistributing it is allowed, and the pieces
are small -- the command-line binary is 1.4 MB and the smallest star database is 100 MB,
less than the offline Gaia catalogue.

**Star databases.** ASTAP matches against its own quad database, which is separate from
the photometric reference catalogue in :mod:`gaiadb` and cannot be substituted for it.
``d05`` (500 stars/deg^2, 100 MB) is the default and covers fields from 0.6 to 10
degrees, which spans both Seestar models -- an S50 sub is 0.72 x 1.28 degrees and an
S30pro 2.24 x 3.99. ``d20`` (400 MB) goes down to 0.3 degrees and is the one to fetch if
a solve is marginal on the S50's narrow axis.

The database is passed to ASTAP explicitly with ``-d``, so it can live in this package's
cache rather than next to the binary or in ``/opt/astap``.
"""

import json
import os
import platform
import stat
import subprocess
import sys
from pathlib import Path

#: Bumped when the mirrored binaries change. The cache directory is keyed on it, so a
#: new version never needs a stale copy invalidated by hand.
VERSION = "2026.07.16"

#: Where the mirrored binaries and databases live. Overridden by ``SEESTAR_ASTAP_URL``,
#: which is how you point at a staging host or a ``file://`` copy without editing code.
#:
#: Mirrored rather than fetched from upstream because SourceForge's download URLs are
#: redirect-heavy and version-stamped, so they rot; MPL-2.0 permits the copy.
BASE_URL = "https://homepage.univie.ac.at/kieran.leschinski/seestar/astap"

#: Filename of the checksum manifest on the mirror.
MANIFEST = "MANIFEST.json"

#: Star database fetched when none is named. See the module docstring.
DEFAULT_DATABASE = "d05"

#: Platform key -> the executable name inside that platform's archive. The key is what
#: :func:`platform_key` resolves to and what the manifest is indexed by.
EXECUTABLES = {
    "win64": "astap_cli.exe",
    "win32": "astap_cli.exe",
    "linux_amd64": "astap_cli",
    "linux_aarch64": "astap_cli",
    "macos_x86_64": "astap_cli",
    "macos_arm64": "astap_cli",
}

#: Names a system-installed ASTAP command-line solver may go by. ``astap_cli`` is what
#: Debian and the AUR install; a manual install is often just renamed ``astap``.
COMMAND_NAMES = ("astap_cli", "astap")


def platform_key():
    """This machine's key into :data:`EXECUTABLES` and the manifest.

    Raises rather than guessing on an unrecognised platform: a wrong binary fails in a
    far more confusing way than a clear refusal.
    """
    machine = platform.machine().lower()
    if sys.platform == "win32":
        return "win64" if machine in ("amd64", "x86_64", "arm64") else "win32"
    if sys.platform == "darwin":
        return "macos_arm64" if machine in ("arm64", "aarch64") else "macos_x86_64"
    if sys.platform.startswith("linux"):
        if machine in ("x86_64", "amd64"):
            return "linux_amd64"
        if machine in ("aarch64", "arm64"):
            return "linux_aarch64"
    raise RuntimeError(
        f"no ASTAP build is mirrored for {sys.platform}/{machine}. Install ASTAP "
        "yourself (https://www.hnsky.org/astap.htm) and it will be found on PATH."
    )


# --- where things live ------------------------------------------------------------------

def cache_dir():
    """Directory the downloaded solver lives in. Downloads nothing.

    Honours ``SEESTAR_ASTAP_DATA`` first, then the platform cache location -- the same
    precedence :func:`examples.data_dir` and :func:`gaiadb.cache_dir` use.
    """
    override = os.environ.get("SEESTAR_ASTAP_DATA")
    if override:
        return Path(override).expanduser()
    base = os.environ.get("XDG_CACHE_HOME") or os.environ.get("LOCALAPPDATA")
    root = Path(base).expanduser() if base else Path.home() / ".cache"
    return root / "seestar-photometry" / f"astap-{VERSION}"


def database_dir():
    """Directory the star databases are unpacked into. Downloads nothing."""
    return cache_dir() / "data"


def downloaded_executable():
    """Path of the downloaded solver, whether or not it exists."""
    try:
        name = EXECUTABLES[platform_key()]
    except RuntimeError:
        return None
    return cache_dir() / name


def has_database(directory=None):
    """Whether a star database is present in the cache. Downloads nothing.

    ASTAP databases are many numbered ``.1476`` files, so presence is judged by finding
    any of them rather than by name.
    """
    directory = Path(directory) if directory is not None else database_dir()
    return directory.is_dir() and any(directory.glob("*.1476"))


def executable(explicit=None):
    """The ASTAP command-line solver to use, or ``None`` if there is none.

    Resolution order, most specific first:

    1. ``explicit`` -- normally ``Project.astap_exe``.
    2. ``$ASTAP_EXE``.
    3. ``astap_cli`` or ``astap`` on ``PATH``. This is what makes a system install work
       with no configuration: Debian and the AUR both ship ``astap-cli``.
    4. A copy fetched by :func:`download`.
    5. The stock Windows install location, which is where the installer puts it and is
       the only platform with one worth guessing.

    Returns a path, not a bare name, so callers can report what they are about to run.
    """
    import shutil

    if explicit:
        found = shutil.which(str(explicit))
        return Path(found) if found else Path(explicit)

    from_env = os.environ.get("ASTAP_EXE")
    if from_env:
        return Path(from_env)

    for name in COMMAND_NAMES:
        found = shutil.which(name)
        if found:
            return Path(found)

    local = downloaded_executable()
    if local is not None and local.exists():
        return local

    stock = Path(r"C:\Program Files\astap\astap_cli.exe")
    return stock if stock.exists() else None


def is_installed():
    """Whether an ASTAP solver can be found at all. Downloads nothing."""
    found = executable()
    return found is not None and Path(found).exists()


# --- fetching -------------------------------------------------------------------------

def _base_url():
    return os.environ.get("SEESTAR_ASTAP_URL", BASE_URL).rstrip("/")


def manifest():
    """The mirror's manifest, fetched fresh. Needs the network."""
    from .gaiadb import _fetch

    directory = cache_dir()
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / MANIFEST
    _fetch(f"{_base_url()}/{MANIFEST}", path)
    return json.loads(path.read_text(encoding="utf-8"))


def _unpack(archive, destination):
    """Extract a zip, flattening it, and return the files written.

    Flattened because the upstream archives are inconsistent about whether they contain
    a top-level directory, and everything in them is wanted in one place anyway.
    """
    import zipfile

    destination.mkdir(parents=True, exist_ok=True)
    written = []
    with zipfile.ZipFile(archive) as zf:
        for member in zf.infolist():
            if member.is_dir():
                continue
            name = Path(member.filename).name
            if not name or name.startswith("."):
                continue
            target = destination / name
            with zf.open(member) as source, open(target, "wb") as fh:
                fh.write(source.read())
            written.append(target)
    return written


def _make_runnable(path):
    """Make a freshly-unpacked binary actually executable.

    Two things stop it. A zip does not carry the Unix execute bit, so it has to be set.
    And on macOS the download is quarantined by Gatekeeper, which refuses to run it at
    all until the attribute is cleared -- ``xattr`` failing is not itself fatal, because
    the verification step is what decides whether this worked.
    """
    if sys.platform == "win32":
        return
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    if sys.platform == "darwin":
        try:
            subprocess.run(["xattr", "-dr", "com.apple.quarantine", str(path)],
                           capture_output=True, timeout=30, check=False)
        except (OSError, subprocess.SubprocessError):
            pass


def verify(path=None):
    """Run the solver once to prove it is usable. Raises with the fix if it is not.

    Worth doing at download time rather than leaving it to the first solve: a blocked or
    non-executable binary otherwise surfaces deep inside a batch run as a per-frame
    failure, once per frame.
    """
    path = Path(path) if path is not None else executable()
    if path is None or not Path(path).exists():
        raise RuntimeError("no ASTAP executable to verify; run astap.download() first")
    try:
        subprocess.run([str(path), "-h"], capture_output=True, timeout=60, check=False)
    except OSError as exc:
        hint = ""
        if sys.platform == "darwin":
            hint = (f"\nmacOS blocks unsigned downloads. Try:\n"
                    f"    xattr -dr com.apple.quarantine {Path(path).parent}\n"
                    "and if that is not enough, approve it once in System Settings > "
                    "Privacy & Security.")
        elif sys.platform != "win32":
            hint = f"\nTry:\n    chmod +x {path}"
        raise RuntimeError(
            f"ASTAP at {path} could not be run: {type(exc).__name__}: {exc}{hint}"
        ) from exc
    return Path(path)


def download(database=DEFAULT_DATABASE, force=False, quiet=False):
    """Fetch the ASTAP solver and a star database, and return the executable's path.

    Idempotent: returns straight away if both are already there. Each piece is
    checksum-verified against the mirror's manifest before it is unpacked, and the
    binary is run once at the end to prove it works.

    Pass ``database=None`` to fetch only the binary -- worth it only if a database is
    already installed somewhere ASTAP looks by itself.
    """
    from .gaiadb import _fetch

    key = platform_key()
    directory = cache_dir()
    target = directory / EXECUTABLES[key]
    want_database = database is not None and not (has_database() and not force)

    if not force and target.exists() and not want_database:
        return target

    meta = manifest()
    if key not in meta.get("platforms", {}):
        raise RuntimeError(
            f"the mirror at {_base_url()} carries no ASTAP build for {key!r} "
            f"(has: {sorted(meta.get('platforms', {}))})"
        )

    if force or not target.exists():
        entry = meta["platforms"][key]
        if not quiet:
            print(f"[seestar-photometry] fetching ASTAP {meta.get('version', '?')} "
                  f"for {key} ({entry['bytes'] / 1e6:.1f} MB)", flush=True)
        archive = directory / entry["file"]
        _fetch(f"{_base_url()}/{entry['file']}", archive, entry.get("sha256"))
        _unpack(archive, directory)
        archive.unlink(missing_ok=True)
        if not target.exists():
            raise RuntimeError(
                f"{entry['file']} did not contain {EXECUTABLES[key]!r}; the mirror's "
                "archive layout may have changed"
            )
        _make_runnable(target)

    if want_database:
        entry = meta.get("databases", {}).get(database)
        if entry is None:
            raise RuntimeError(
                f"the mirror carries no star database {database!r} "
                f"(has: {sorted(meta.get('databases', {}))})"
            )
        if not quiet:
            print(f"[seestar-photometry] fetching the {database} star database "
                  f"({entry['bytes'] / 1e6:.0f} MB)", flush=True)
        archive = directory / entry["file"]
        _fetch(f"{_base_url()}/{entry['file']}", archive, entry.get("sha256"))
        _unpack(archive, database_dir())
        archive.unlink(missing_ok=True)

    verify(target)
    if not quiet:
        print(f"[seestar-photometry] ASTAP ready: {target}", flush=True)
    return target
