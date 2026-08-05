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
| `dev` | the above, plus `pytest` | Running the test suite |

```bash
pip install "seestar-photometry[plot,stack]"
```

`matplotlib` is an extra rather than a dependency because the measurement path is imported
by batch jobs on headless machines; it is imported lazily, so `import seestar_photometry`
never needs it. In practice you almost certainly want it — the diagnostics are how you
check a reduction.

## A plate solver

The on-board Seestar WCS is not accurate enough for photometry (see
[astrometry-and-gaia](astrometry-and-gaia.md)), so frames must be re-solved. Two options:

**ASTAP** — the default. Local, offline, about a second per frame, no key, no rate limit.
[Download](https://www.hnsky.org/astap.htm) it plus a star database (`D50` is plenty), then
point the package at the binary if it is not in the default location:

```python
Project(..., astap_exe=r"C:\Program Files\astap\astap_cli.exe")
```

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
