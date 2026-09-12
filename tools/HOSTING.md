# Hosting the downloadable data

Three things are fetched at runtime rather than shipped in the wheel. Two of them live on
the CrowdSky webspace; the third is a GitHub release asset.

| What | Size | Where | Built by |
|---|---|---|---|
| Example data | 18 MB | GitHub release `example-data-v1` | `tools/build_example_data.py` |
| Gaia catalogue | 6.35 GB, all sky | `crowdsky.univie.ac.at/seestar_assets/gaia-seestar-dr3-v2/` | `build_gaia_allsky.py` + `build_gaia_s50.py` |
| ASTAP mirror | 106 MB | `crowdsky.univie.ac.at/seestar_assets/astap/` | `tools/build_astap_mirror.py` |

The CrowdSky webspace is mounted over the VPN at `Z:\crowdskyo92\html\seestar_assets`,
which maps to `https://crowdsky.univie.ac.at/seestar_assets/`. Deploying is a file copy;
there is no rsync-over-ssh step. **50 GB total**, so the 6.35 GB catalogue is the one asset
whose size is worth thinking about before adding another.

> **The share refuses attribute and timestamp writes.** A plain `robocopy` fails on every
> directory with `ERROR 50 ... request not supported`, having copied nothing. Copy the data
> and nothing else -- `/COPY:D /DCOPY:D` -- and it works. The same goes for anything else
> that tries to preserve metadata; `cp -p` and `rsync -a` will fail the same way.
> Note also that robocopy exit codes **below 8 are success** (1 means "files were copied"),
> so a wrapper that tests for zero will report a good copy as a failure.

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
  gaia-seestar-dr3-v2/
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

**Currently deployed: the full sky, as `v2`.** 12288 parts, 148,667,792 rows, 6.35 GB. It
replaces nothing -- until August 2026 the server carried no Gaia catalogue at all, only the
ASTAP mirror, so every client was silently falling back to TAP. A region-only build (MW Cam,
26 parts, 2.5 MB) existed locally and was never uploaded.

The manifest is re-fetched on every call, so clients pick up a wider or corrected dataset
without being told anything.

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

**Eleven of the schema's columns come out populated**, and they are the ones the
package reads: `source_id`, `ra`, `dec`, `pmra`, `pmdec`, `phot_g_mean_mag`, `v_jkc_mag`,
`b_jkc_mag`, `r_jkc_mag`, `v_jkc_flag`, `c_star`. Absent, because they live only in
`gaia_source`: `bp_rp`, `phot_variable_flag`, `ruwe`, `teff_gspphot`,
`ipd_frac_multi_peak`, `non_single_star`, `duplicated_source`, `in_galaxy_candidates`,
`has_epoch_photometry`. The six `*_s50_mag` columns come from the separate pass below.
Everything unfilled is written as an entirely-null column, so it reads back masked and the
schema is the same either way; `columns_present` in the manifest says which are real.

The one that costs something is `phot_variable_flag`: `fit_zeropoint` and
`select_comparisons` both use it to drop known variables, and without it they drop none.
Both guard on the column being present, so nothing breaks. It could be recovered
cheaply -- DR3 has ~10.5M variable sources and only their `source_id` is needed -- which is
the obvious next improvement. `bp_rp` is not a loss: the comparison colour cut falls back
to B-R, which this build does have.

Bump `gaiadb.DATA_VERSION` when the schema or the row selection changes; bump
`DATA_RELEASE` for a new Gaia release. The directory name carries both, so two versions can
sit side by side on the server during a migration.

### Synthetic S50 magnitudes (`v2`)

`tools/build_gaia_s50.py` adds `r_s50_mag` / `g_s50_mag` / `b_s50_mag` and their errors by
integrating Gaia DR3 XP spectra through the Seestar's fitted passbands
(`tools/bayer_passbands.ecsv`, copied from `photometry-adventures-2`, B red tail ×1.15).
Calibrating a Bayer plane against its own band removes the zero point's colour term instead
of fitting it — measured on 57 fifteen-minute 20 Cep stacks, G +0.054 → −0.002, R −0.132 →
−0.002, B +0.548 → ≈0 mag/mag, and the B-plane calibrator rms fell from 347 to 104 mmag.
This is the `v1` → `v2` schema change, so **`v1` must stay on the server until `v2` is
uploaded**: `DATA_VERSION` is in the URL and clients follow it.

It adds columns to an existing dataset rather than rebuilding one, so `--src` is the `v1`
tree and `--out` the `v2` tree:

```bash
uv run python tools/build_gaia_s50.py --stage list --work D:/tmp/xp
uv run python tools/build_gaia_s50.py --validate --work D:/tmp/xp
uv run python tools/build_gaia_s50.py --stage all --work D:/tmp/xp \
    --src D:/tmp/gaia_allsky --out D:/tmp/gaia_s50 --workers 8
```

**The spectra stop at G = 15.0.** ESA's bulk `xp_sampled_mean_spectrum` product is the
bright subset of the XP population, not all of it. Measured over eleven fully covered
HEALPix parts: coverage of the V < 17.5 catalogue is **100% below G = 15.0 and exactly 0%
above it**, which is 23.03% of all rows — 37% in a high-latitude part, 17% in the Galactic
plane, where the same G cut takes a smaller share. The other three quarters get masked
columns and keep calibrating against `v_jkc_mag`. Full depth would mean
`xp_continuous_mean_spectrum` instead: **3.6 TB** against 114 GB, plus a GaiaXPy sampling
pass per source rather than a dot product. Worth revisiting only if the faint end turns out
to matter.

**The download is the entire job, and it needs concurrency to go anywhere.** ESA throttles
a *connection* to ~4 MB/s, not a client, so one stream fetches 114.4 GB in eight hours and
sixteen fetch it in forty minutes. Measured over a wired link from Vienna: 4.4 MB/s on one
stream, 50 on 16, then flat — 46 at 24, 50 at 32, 44 at 64. Sixteen is the default because
that is where it plateaus.

Check the link before blaming the CDN. The same probe over WiFi gave 3.3 MB/s on one stream
and 5.1 on 48, which reads exactly like a server-side per-client cap and is not one — it was
the local link saturating at ~40 Mbit/s. `--download-workers 1` against `--download-workers
16` tells the two apart in a minute.

The synthesis stage does not care where it runs: **71 minutes for all 12288 parts on eight
workers**, at ~1.2 GB resident each. So the whole job is about two hours on any machine with
a real uplink, and moving it to a bigger host buys nothing.

Measured on the full build: 148,667,792 rows — the same count as `v1`, since the stage adds
columns and never filters — 34,238,688 of them with S50 magnitudes, in 12288 parts totalling
**6.35 GB** against `v1`'s 5.65 GB. 43.7 bytes a row. Six float32 columns that are three
quarters null cost 0.70 GB, which is what byte-stream-split plus zstd does with them.

Both stages are resumable. `fetch` skips any file already on disk at the size the listing
gives; `synth` skips any part already written, so a re-run after a kill costs seconds.
Keep the 114 GB afterwards if there is room — it is reusable for any re-synthesis, and
re-fetching it is another day.

Validation before adopting a build:

```bash
uv run python tools/build_gaia_s50.py --validate --work D:/tmp/xp
```

recomputes the 20 Cep field and compares against
`photometry-adventures-2/paper/data/gaia_20Cep_bayer.ecsv`. Identical inputs and formula,
so it should agree exactly; it reports the worst disagreement and exits non-zero above
0.0001 mag. Pass a path to `--validate` to check against a different reference table.

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
