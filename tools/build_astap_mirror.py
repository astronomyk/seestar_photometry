#!/usr/bin/env python
"""Assemble the ASTAP mirror that :func:`astap.download` fetches from.

ASTAP is `MPL-2.0 <https://github.com/han-k59/astap>`_, so mirroring it is allowed. It is
mirrored rather than linked because SourceForge's download URLs are redirect-heavy and
version-stamped: they rot, and a rotted URL breaks the feature for every user at once.

This takes files you have already downloaded rather than fetching them itself. Upstream
publishes the command-line builds on hnsky.org and the star databases on SourceForge,
under names and layouts that change between releases -- pinning that here would be one
more thing to rot. Downloading them by hand once per ASTAP release is the honest cost.

Get the pieces from:

* https://www.hnsky.org/astap.htm -- ``astap_command-line_version_*.zip``, one per
  platform, 1.4 MB each
* https://sourceforge.net/projects/astap-program/files/star_databases/ -- ``d05`` (100
  MB) is the default and covers 0.6-10 degree fields, which spans both Seestar models;
  ``d20`` (400 MB) reaches down to 0.3 degrees

Then::

    uv run python tools/build_astap_mirror.py --out mirror \\
        --version 2026.07.16 \\
        --binary win64=astap_command-line_version_win64.zip \\
        --binary linux_amd64=astap_command-line_version_Linux_amd64.zip \\
        --binary macos_x86_64=astap_command-line_version_macOS_x86_64.zip \\
        --binary macos_arm64=astap_command-line_version_macOS_M1.zip \\
        --database d05=d05_star_database.zip

and upload the contents of ``mirror/`` to the path in ``astap.BASE_URL``. The whole
deployment, including where the files go and how to verify one, is in ``tools/HOSTING.md``.

Every archive is checked for the executable it is supposed to contain before it is
accepted, so a mislabelled platform fails here rather than on a user's machine.

Maintenance script, not shipped in the wheel.
"""

import argparse
import json
import shutil
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from seestar_photometry import astap  # noqa: E402


def _pair(text):
    """``key=path`` as a two-tuple."""
    if "=" not in text:
        raise argparse.ArgumentTypeError(f"expected key=path, got {text!r}")
    key, path = text.split("=", 1)
    return key, Path(path)


def _sha256(path, chunk=1 << 20):
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def _check_binary(archive, key):
    """Refuse an archive that does not hold the executable that platform needs."""
    wanted = astap.EXECUTABLES[key]
    with zipfile.ZipFile(archive) as zf:
        names = {Path(n).name for n in zf.namelist()}
    if wanted not in names:
        raise SystemExit(
            f"{archive} does not contain {wanted!r} (it has: {sorted(names)[:8]}).\n"
            f"Is it really the {key} build?"
        )


def _check_database(archive):
    """Refuse an archive with no star-database files in it."""
    with zipfile.ZipFile(archive) as zf:
        if not any(n.endswith(".1476") for n in zf.namelist()):
            raise SystemExit(
                f"{archive} contains no *.1476 files, so it is not an ASTAP star "
                "database. The Windows downloads are self-extracting .exe files -- "
                "unpack one and re-zip the .1476 files, or take the Linux archive."
            )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True, help="mirror directory to build")
    ap.add_argument("--version", default=astap.VERSION,
                    help="ASTAP version being mirrored; must match astap.VERSION")
    ap.add_argument("--binary", type=_pair, action="append", default=[],
                    metavar="KEY=ZIP", help=f"one of {sorted(astap.EXECUTABLES)}")
    ap.add_argument("--database", type=_pair, action="append", default=[],
                    metavar="NAME=ZIP", help="e.g. d05=d05_star_database.zip")
    args = ap.parse_args()

    if not args.binary:
        raise SystemExit("at least one --binary is required")
    if args.version != astap.VERSION:
        print(f"warning: --version {args.version} != astap.VERSION {astap.VERSION}; "
              "bump the module constant or the cache will not be re-keyed", flush=True)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta = {"version": args.version, "platforms": {}, "databases": {}}

    for key, source in args.binary:
        if key not in astap.EXECUTABLES:
            raise SystemExit(f"unknown platform {key!r}; "
                             f"expected one of {sorted(astap.EXECUTABLES)}")
        _check_binary(source, key)
        target = out / f"astap_cli_{key}.zip"
        shutil.copy(source, target)
        meta["platforms"][key] = {"file": target.name, "sha256": _sha256(target),
                                  "bytes": target.stat().st_size}
        print(f"  {key:14s} {target.name:32s} {target.stat().st_size/1e6:7.1f} MB")

    for name, source in args.database:
        _check_database(source)
        target = out / f"{name}.zip"
        shutil.copy(source, target)
        meta["databases"][name] = {"file": target.name, "sha256": _sha256(target),
                                   "bytes": target.stat().st_size}
        print(f"  {name:14s} {target.name:32s} {target.stat().st_size/1e6:7.1f} MB")

    (out / astap.MANIFEST).write_text(json.dumps(meta, indent=1), encoding="utf-8")
    missing = sorted(set(astap.EXECUTABLES) - set(meta["platforms"]))
    if missing:
        print(f"\nnote: no build mirrored for {missing} -- users on those platforms "
              "will be told to install ASTAP themselves")
    if astap.DEFAULT_DATABASE not in meta["databases"]:
        print(f"note: the default database {astap.DEFAULT_DATABASE!r} is not mirrored; "
              "astap.download() will fail unless callers name another")
    print(f"\nwrote {out / astap.MANIFEST}")
    print(f"upload the contents of {out} to {astap.BASE_URL}")


if __name__ == "__main__":
    main()
