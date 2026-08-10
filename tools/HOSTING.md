# Hosting the downloadable data

Three things are fetched at runtime rather than shipped in the wheel. Two of them live on
the CrowdSky webspace; the third is a GitHub release asset.

| What | Size | Where | Built by |
|---|---|---|---|
| Example data | 18 MB | GitHub release `example-data-v1` | `tools/build_example_data.py` |
| Gaia catalogue | ~10 GB full sky | `crowdsky.univie.ac.at/seestar_assets/gaia-seestar-dr3-v1/` | `tools/build_gaia_catalogue.py` |
| ASTAP mirror | 106 MB | `crowdsky.univie.ac.at/seestar_assets/astap/` | `tools/build_astap_mirror.py` |

The CrowdSky webspace is mounted over the VPN at `Z:\crowdskyo92\html\seestar_assets`,
which maps to `https://crowdsky.univie.ac.at/seestar_assets/`. Deploying is a file copy;
there is no rsync-over-ssh step.

The URLs are `gaiadb.BASE_URL` and `astap.BASE_URL`. Both honour an environment override
(`SEESTAR_GAIA_URL`, `SEESTAR_ASTAP_URL`), which is how to test against a staging copy — or
a local directory over `file://` — without editing code.

Everything is served as plain static files over HTTPS. No range requests, no directory
listing, no server-side anything: each download is a whole file fetched by name from a
manifest, so any web server will do.

## Layout on the server

```
Z:\crowdskyo92\html\seestar_assets\        ->  https://crowdsky.univie.ac.at/seestar_assets/
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
2. Get a star database. `d05` is the default: 1476 files, 109 MB raw, 102 MB zipped.
   **The Windows downloads from
   <https://sourceforge.net/projects/astap-program/files/star_databases/> are
   self-extracting `.exe` files**, so the easiest source is an existing ASTAP install --
   zip its `d05_*.1476` files directly. `build_astap_mirror.py` refuses an archive with no
   `.1476` files in it rather than letting a broken one reach the server.

   **`d05` alone is enough for both Seestar models.** Worth stating because an install
   holding both `d05` and `d20` will have ASTAP silently preferring `d20`, so local
   testing can flatter a `d05`-only mirror. Verified against the mirror with no other
   database present: 6 of 6 real MW Cam raw subs solved, three per model, in 0.2-0.4 s.
   `d20` is 423 MB and buys nothing here.
3. Build and upload:

   ```bash
   uv run python tools/build_astap_mirror.py --out mirror/astap \
       --binary win64=astap_command-line_version_win64.zip \
       --binary linux_amd64=astap_command-line_version_Linux_amd64.zip \
       --binary macos_x86_64=astap_command-line_version_macOS_x86_64.zip \
       --binary macos_arm64=astap_command-line_version_macOS_M1.zip \
       --database d05=d05.zip

   cp -r mirror/astap "Z:/crowdskyo92/html/seestar_assets/"
   ```

   Measured: 106 MB copied over the VPN in 37 s.

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

Copy the result into `Z:\crowdskyo92\html\seestar_assets\` and it is live.

**Fetch the radius the project needs, not the frame's.** `Project.catalogue_radius_deg`
(2.221 degrees at the default `catalogue_half_deg=1.5`) reaches the corners of the cached
box, and the cached box covers a whole night's dithering. Asking for 1.5 degrees pulls 11
HEALPix parts where the project wants 14, and the only symptom is
`catalogue_backend_used()` quietly answering `"tap"`. `Project.fetch_catalogue()` gets this
right, so prefer it:

```python
proj.fetch_catalogue()          # not gaiadb.download(radius_deg=<fov>)
```

**Currently deployed: MW Cam only.** A 4-degree region, 47995 rows in 26 parts, 2.5 MB --
enough to exercise the offline path end to end on the main science field. Everywhere else
`covers()` returns False and `catalogue_backend="auto"` falls back to TAP, which is the
designed behaviour and not an error. Replace it with a wider build when there is one; the
manifest is re-fetched on every call, so clients pick up new regions without being told.

### All sky

`tools/build_gaia_allsky.py`. The obvious route does not work, and the reasons are worth
keeping because both are counter-intuitive.

**Bulk-downloading everything is not viable.** The catalogue needs two Gaia tables.
`synthetic_photometry_gspc` has the synthetic V and the colour and bulk-downloads as 44
files totalling ~42 GB. `gaia_source` has the positions and proper motions and
bulk-downloads as 3387 files totalling **~790 GB** -- 152 columns for 1.8 billion rows, to
keep 6 columns for the 220 million that have synthetic photometry. And neither table alone
is enough: **GSPC carries no `ra`/`dec`**, only `source_id`.

**Chunking the join over TAP does not work either.** `source_id BETWEEN` has no spatial
index behind it; the archive answers HTTP 500 after exactly 182 seconds, measured three
times, on queries as simple as a row count. Cone-tiling works but means ~525 queries.

So the build **joins new photometry onto positions already in hand**. Given per-HEALPix
tiles carrying `source_id`/`ra`/`dec`/`pmra`/`pmdec` -- `D:	mp\gaia_V
side32` is such a
set, 12288 tiles at nside 32, which is exactly `gaiadb`'s partitioning -- only the 42 GB of
GSPC has to be fetched:

```bash
uv run python tools/build_gaia_allsky.py     --tiles D:/tmp/gaia_V/nside32 --work D:/tmp/gaia_build --out D:/tmp/gaia_allsky
```

Two passes, both resumable: a routed GSPC file leaves a marker and pass 2 skips parts
already written, so it survives being killed or the laptop sleeping. Downloads run four at
a time because the CDN gives ~2 MB/s on one stream and ~5 MB/s over four; routing stays
single-threaded, which is why the bucket appends need no locking. Reckon on ~2.5 h of
download, ~12 min of parsing and well under an hour to join.

**Eleven of the twenty schema columns come out populated**, and they are the ones the
package reads: `source_id`, `ra`, `dec`, `pmra`, `pmdec`, `phot_g_mean_mag`, `v_jkc_mag`,
`b_jkc_mag`, `r_jkc_mag`, `v_jkc_flag`, `c_star`. Absent, because they live only in
`gaia_source`: `bp_rp`, `phot_variable_flag`, `ruwe`, `teff_gspphot`,
`ipd_frac_multi_peak`, `non_single_star`, `duplicated_source`, `in_galaxy_candidates`,
`has_epoch_photometry`. They are written as entirely-null columns, so they read back masked
and the schema is the same either way; `columns_present` in the manifest says which are
real.

The one that costs something is `phot_variable_flag`: `fit_zeropoint` and
`select_comparisons` both use it to drop known variables, and without it they drop none.
Both guard on the column being present, so nothing breaks. It could be recovered
cheaply -- DR3 has ~10.5M variable sources and only their `source_id` is needed -- which is
the obvious next improvement. `bp_rp` is not a loss: the comparison colour cut falls back
to B-R, which this build does have.

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
