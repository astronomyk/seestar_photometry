# Installation

```bash
pip install seestar-photometry
```

Python 3.11 or newer. The core install pulls `astropy`, `astroquery`, `numpy`, `scipy` and
`sep` — about 100 kB of package on top of those.

## The example data

Every example in these docs runs on real Seestar frames, but they are **not** in the wheel:
18 MB of FITS in a package everyone installs would be rude to PyPI and to anyone who only
wants the photometry. They are a GitHub release asset, fetched once on first use and cached:

```python
from seestar_photometry import examples

examples.download()          # ~18 MB, once
frame = examples.stack()     # every accessor fetches implicitly if needed
```

The cache location, in order of precedence:

| | |
|---|---|
| `SEESTAR_PHOTOMETRY_DATA` | an explicit directory — set this and nothing is downloaded |
| `XDG_CACHE_HOME` / `LOCALAPPDATA` | the platform cache location |
| otherwise | `~/.cache/seestar-photometry` |

**Offline or air-gapped?** Fetch
[`example-data-v1.tar.gz`](https://github.com/astronomyk/seestar_photometry/releases/download/example-data-v1/example-data-v1.tar.gz)
elsewhere, unpack it, and point `SEESTAR_PHOTOMETRY_DATA` at the directory. The archive is
checksum-verified on download, so a truncated file fails loudly rather than surfacing later
as inexplicable photometry.

In CI, cache that directory keyed on the dataset version — see `.github/workflows/tests.yml`.

## Extras

| Extra | Adds | You want it for |
|---|---|---|
| `plot` | `matplotlib` | Any diagnostic figure — `plots`, `report` |
| `stack` | `astroalign`, `scikit-image` | Stacking raw sub-exposures |
| `catalog` | `pyarrow`, `astropy-healpix` | The offline Gaia catalogue and `solver="local"` |
| `dev` | the above, plus `pytest` | Running the test suite |

```bash
pip install "seestar-photometry[plot,stack,catalog]"
```

`matplotlib` is an extra rather than a dependency because the measurement path is imported
by batch jobs on headless machines; it is imported lazily, so `import seestar_photometry`
never needs it. In practice you almost certainly want it — the diagnostics are how you
check a reduction.

## The offline Gaia catalogue

Optional, and much larger than the example data. It removes the Gaia TAP query — the
least reliable step in the pipeline — and it is what `solver="local"` pairs against.

```bash
pip install "seestar-photometry[catalog]"
```

```python
from seestar_photometry import gaiadb

gaiadb.download(center=(186.68, 81.47), radius_deg=5)   # just this patch of sky
gaiadb.download()                                        # or the whole thing
```

Fetching per region is the point of the layout: the catalogue is partitioned into 12288
HEALPix tiles, so a few observing fields cost tens of MB rather than the full multi-GB
set. Its own cache location, so it can live on a different disk from the example data:

| | |
|---|---|
| `SEESTAR_GAIA_DATA` | an explicit directory — set this and nothing is downloaded |
| `XDG_CACHE_HOME` / `LOCALAPPDATA` | the platform cache location |
| otherwise | `~/.cache/seestar-photometry` |

Nothing has to be configured to use it. `Project` defaults to
`catalogue_backend="auto"`, which uses the local copy when it covers the field and falls
back to a TAP query when it does not — so installing the download is the only step.
Force one or the other with `catalogue_backend="local"` (raises rather than reaching the
network, which is what you want on a machine that must not) or `"tap"`.

It also carries proper motions, which the TAP query does not fetch at all. Gaia
positions are J2016.0, and by 2026 a 100 mas/yr star has moved 1 arcsec — half the
default match tolerance:

```python
Project(..., epoch=2026.4)      # propagate to the observing season
```

To build a catalogue for one field yourself instead of downloading, from a source
checkout:

```bash
uv run python tools/build_gaia_catalogue.py --ra 186.6821 --dec 81.474 --radius 3
```

## A plate solver

The on-board Seestar WCS is not accurate enough for photometry (see
[astrometry-and-gaia](astrometry-and-gaia.md)), so frames must be re-solved. Three
options:

**The local solver** — no binary, no network, no index files. It pairs detections against
the reference catalogue the pipeline already needs, so the only requirement is having that
catalogue: either the offline download above, or a cached TAP one from a previous run.

```python
Project(..., solver="local")
```

Because it is anchored rather than blind it needs the frame's header pointing, which every
Seestar writes. For a frame with no pointing at all, use ASTAP.

Accuracy is comparable to ASTAP — on real MW Cam data the two agree to 0.5–2 arcsec, well
under a pixel — but it is **several times slower** (1–7 s per frame against 0.2–0.4 s), and
on an S30pro a minority of raw subs fail to solve at all. They fail loudly, so the batch
runner records them and moves on. Pick this one when you would rather not install anything;
pick ASTAP when you are solving thousands of frames.

**ASTAP** — the default, and the fastest. Local, offline, no key, no rate limit.

You do not have to configure it. The package looks for `astap_cli` or `astap` on `PATH`
before anything else, so a system install just works:

```bash
apt install astap-cli          # Debian/Ubuntu; also in the AUR
```

If you would rather not install a system package, it can fetch its own copy — a 1.4 MB
binary and a 100 MB star database, into the same cache the other downloads use:

```python
from seestar_photometry import astap

astap.download()               # once
astap.executable()             # -> what solve_astap will run
```

The star database is ASTAP's own quad index and is **not** the same thing as the Gaia
reference catalogue above; you need both if you use this solver. `d05` is the default and
covers 0.6–10° fields, which spans both Seestar models. Pass `database="d20"` (400 MB) if a
solve is marginal on the S50's narrow axis.

Resolution order, if you want to override it: `Project(astap_exe=...)` → `$ASTAP_EXE` →
`PATH` → the fetched copy → `C:\Program Files\astap\astap_cli.exe`.

On macOS the downloaded binary is unsigned, so Gatekeeper quarantines it. `astap.download()`
clears the quarantine attribute and then runs the binary once to check; if it is still
blocked you get an error naming the command to run. This path is untested on macOS — a
system install via [the official download](https://www.hnsky.org/astap.htm) avoids it.

**astrometry.net** — the fallback, for frames ASTAP cannot solve. Set an API key:

```bash
export ASTROMETRY_KEY=xxxxxxxxxxxx
```

```python
Project(..., solver="nova")
```

The key is read from the environment and never stored. Note the service accepts only a
handful of concurrent jobs per account, so `solve_all` caps itself at four threads for this
solver.

Neither is needed to follow these docs: the example frames come with solved `.wcs` sidecars.

## Verifying the install

```python
from seestar_photometry import examples

frame = examples.stack()      # fetches the example data on first call
print(frame.layout, frame.shape, frame.model)
# cube (1000, 1000) S50
```

To run the test suite from a source checkout:

```bash
git clone https://github.com/astronomyk/seestar_photometry
cd seestar_photometry
uv sync --extra dev
uv run pytest
```

The suite downloads the example data on first run (and skips the real-data module entirely if
that is not possible, so an offline checkout still tests the algorithms). It combines
synthetic frames with injected Gaussian PSFs of
known flux — so it asserts *recovery* of a known zero point, colour term, aperture and
period — with tests on the real example frames, which is where instrument surprises show up.
