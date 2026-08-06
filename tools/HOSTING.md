# Hosting the downloadable data

Three things are fetched at runtime rather than shipped in the wheel. Two of them live on
the CrowdSky webspace; the third is a GitHub release asset.

| What | Size | Where | Built by |
|---|---|---|---|
| Example data | 18 MB | GitHub release `example-data-v1` | `tools/build_example_data.py` |
| Gaia catalogue | ~10 GB full sky | `crowdsky.univie.ac.at/seestar/gaia-seestar-dr3-v1/` | `tools/build_gaia_catalogue.py` |
| ASTAP mirror | ~102 MB | `crowdsky.univie.ac.at/seestar/astap/` | `tools/build_astap_mirror.py` |

The URLs are `gaiadb.BASE_URL` and `astap.BASE_URL`. Both honour an environment override
(`SEESTAR_GAIA_URL`, `SEESTAR_ASTAP_URL`), which is how to test against a staging copy — or
a local directory over `file://` — without editing code.

Everything is served as plain static files over HTTPS. No range requests, no directory
listing, no server-side anything: each download is a whole file fetched by name from a
manifest, so any web server will do.

## Layout on the server

```
seestar/
  astap/
    MANIFEST.json
    astap_cli_win64.zip            1.4 MB
    astap_cli_linux_amd64.zip
    astap_cli_macos_x86_64.zip
    astap_cli_macos_arm64.zip
    d05.zip                        100 MB
  gaia-seestar-dr3-v1/
    MANIFEST.json
    hpx5=00000/part.parquet
    ...                            12288 of these
```

Both manifests carry a SHA-256 and a byte count per file, and both downloaders verify
before moving anything into place. A truncated transfer therefore fails loudly instead of
surfacing months later as a mysteriously empty patch of sky.

## ASTAP mirror

ASTAP is [MPL-2.0](https://github.com/han-k59/astap), so mirroring the binaries is allowed.
It is mirrored rather than hot-linked because SourceForge's download URLs are version-stamped
and rot.

1. Download the command-line builds from <https://www.hnsky.org/astap.htm> — one
   `astap_command-line_version_*.zip` per platform, 1.4 MB each.
2. Download a star database from
   <https://sourceforge.net/projects/astap-program/files/star_databases/>. `d05` (100 MB) is
   the default and covers 0.6-10 degree fields, which spans both Seestar models.
   **The Windows downloads are self-extracting `.exe` files**: unpack one and re-zip the
   `.1476` files, or take the Linux archive. `build_astap_mirror.py` refuses an archive with
   no `.1476` files in it rather than letting a broken one reach the server.
3. Build and upload:

   ```bash
   uv run python tools/build_astap_mirror.py --out mirror/astap \
       --binary win64=astap_command-line_version_win64.zip \
       --binary linux_amd64=astap_command-line_version_Linux_amd64.zip \
       --binary macos_x86_64=astap_command-line_version_macOS_x86_64.zip \
       --binary macos_arm64=astap_command-line_version_macOS_M1.zip \
       --database d05=d05.zip

   rsync -av mirror/astap/ crowdsky.univie.ac.at:/var/www/seestar/astap/
   ```

Platforms you do not mirror are not broken — users there are told to install ASTAP
themselves, which on Debian and Arch is one `apt`/`pacman` away.

When ASTAP releases a new version, bump `astap.VERSION` as well. The cache directory is
keyed on it, so a bump is what makes clients pick up the new binary instead of the one they
already have.

## Gaia catalogue

`build_gaia_catalogue.py` queries one region at a time and merges regions into a single
directory, so a dataset is built by running it repeatedly:

```bash
uv run python tools/build_gaia_catalogue.py --ra 186.68 --dec 81.47 --radius 3 --out build/gaia
```

That is the whole story for a per-field or per-project dataset, and it is enough to use the
offline path today.

**A full-sky build is a different job and is not automated.** At V < 17.5 the
synthetic-photometry catalogue is ~195 million rows. Two ways to get there:

- **Tile the sky with TAP cones.** ~525 five-degree cones, merged into one output directory.
  It works with the tool as written and needs no new code, but it is a long unattended run
  and Gaia TAP truncates large async results, so expect to re-run failed tiles.
- **Bulk download and join offline.** The archive publishes `gaia_source` and
  `synthetic_photometry_gspc` as per-HEALPix files under `cdn.gea.esac.esa.int`. Far faster
  and far more reliable, but it needs a fetch-and-join step that does not exist yet. Check
  the current directory layout before scripting it — that CDN is a JavaScript file browser
  and its paths have moved between releases.

Either way the writing half is already shared: `gaiadb.write_dataset` takes whatever table it
is handed, computes the HEALPix partition from `source_id`, and emits the manifest.

Expect roughly 10 GB. Measured on a real 4-degree region: 47995 rows in 26 parts at 51.4
bytes/row, though per-part overhead is inflated at that size and the full build should come
in lower.

Bump `gaiadb.DATA_VERSION` when the schema or the row selection changes; bump
`DATA_RELEASE` for a new Gaia release. The directory name carries both, so two versions can
sit side by side on the server during a migration.

## Verifying a deployment

From a machine that has never fetched anything:

```bash
python -c "
from seestar_photometry import astap
astap.download()          # binary + d05, checksum-verified, then run once
print(astap.executable())
"

python -c "
from seestar_photometry import gaiadb
gaiadb.download(center=(186.68, 81.47), radius_deg=3)
t = gaiadb.cone((186.68, 81.47), 1.5)
print(len(t), 'sources', t.colnames[:4])
"
```

`astap.download()` runs the binary once before returning, so a successful call means the
mirror is serving a working solver and not just bytes. To rehearse without touching the real
server, point the override at a local build:

```bash
SEESTAR_ASTAP_URL="file:///path/to/mirror/astap" python -c "..."
```

## macOS

The mirrored macOS binaries are unsigned, so Gatekeeper quarantines them on download.
`astap.download()` clears the quarantine attribute and then runs the binary to check, raising
with the exact command to run if it is still blocked. **This path has never been tested on a
real Mac.** Worth one run on Apple silicon before telling anyone it works; if Gatekeeper
proves stubborn, the honest fix is to drop the macOS entries from the mirror and let those
users install ASTAP themselves.
