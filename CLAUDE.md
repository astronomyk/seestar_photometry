# seestar-photometry — conventions

Read this before changing code. It is the *shape* of the package; the rationale for every
numerical choice lives in `docs/`, and the API reference is the docstrings.

## Names

Distribution `seestar-photometry`, import package `seestar_photometry`, `src/` layout,
hatchling backend. Python import names can't contain hyphens.

## Commands

`uv` only. Never `uv run --active`. Point PyCharm at `.venv`; don't let it make its own.

```bash
uv sync --extra dev          # matplotlib, pytest, astroalign, scikit-image, build, twine
uv run pytest                # whole suite, offline, ~seconds
uv run pytest tests/test_photometry.py -k aperture
uv run python examples/mwcam_lightcurve.py
uv run sphinx-build -b html -W docs docs/_build/html   # -W: CI builds docs warning-free
```

CI (`.github/workflows`) runs the suite on ubuntu + windows × py3.11/3.13, and builds the
docs with `-W`. There is no linter or formatter configured — match the surrounding style.

## Hard conventions

- **`BANDS = ("R", "G", "B")` is the canonical axis-0 order.** Every frame is a
  `(3, ny, nx)` float32 array in that order. Every per-band array (`Extraction.rms`,
  `.aperture`, `.fwhm`) is shape `(3,)` in that order. **Green is index 1 and is the
  science band** (≈ Johnson V through the Seestar's IRCUT filter). Only green is
  calibrated onto reference V; R and B stay instrumental, which is all the colour
  diagnostics need.

- **Image shape comes from the data, never the header.** `frame.shape` is the accessor.
  The CrowdSky multi-extension layout has an *empty* primary HDU with no `NAXIS1`/`NAXIS2`,
  so any code reading the shape from a header is broken on half the datasets.

- **The on-board WCS is off by ~1 arcmin and unusable for photometry.** Every frame is
  re-solved and cached as a `.wcs` sidecar. Don't "optimise" by reading the header WCS.
  `astrometry.lift` is the deliberate exception, for CrowdSky frames solved server-side.

- **`SeestarFrame` and `Project` are thin data containers, not god objects.** All
  analysis is module-level *functions* taking the object as their first argument
  (`photometry.extract_sources(frame, ...)`). Methods that briefly lived on the frame
  class were moved out on purpose. Configuration is Python on `Project`, not files —
  retargeting means constructing a different `Project`, never editing code.

- **Optional dependencies are imported lazily**, inside the function that needs
  them — never at module top. `matplotlib` in particular: `import seestar_photometry` must
  work in a core install. `tests/test_imports.py` asserts `matplotlib` *and* `astroquery`
  are absent from `sys.modules` after import (in a fresh interpreter), so this can't rot.
  The same applies to `astroalign` / `scikit-image` (the `stack` extra).

- **The standard photometry path and the cross-checks live in separate modules.**
  `photometry` is what the science uses; `kron` (Kron/AUTO) is kept only for
  star/galaxy separation and total-flux cross-checks. Do not promote it.

- **Derived products never go into the data tree or the FITS headers.** Tables and
  figures go to `Project.work_dir`. The one exception is the per-frame `.wcs`
  sidecar, which lives next to its FITS file because it is expensive to recompute
  and shared between projects.

- **No defensive error handling in the measurement layer.** If a frame is broken,
  let it raise; the batch runner in `pipeline` catches per-frame exceptions, records
  a status, and keeps going. Don't scatter `try/except` through the science code.

## Module layers

Each module imports only from layers above it.

| Layer | Modules |
|---|---|
| data | `frames`, `debayer` |
| build | `stacking` (raw subs → a frame; imports only `frames`/`debayer`) |
| measure | `photometry`, `astrometry`, `catalogs` |
| calibrate | `calibration`, `quality` |
| science | `lightcurves`, `contamination`, `depth` |
| orchestrate | `project`, `pipeline` |
| inspect | `plots`, `report`, `_style` |
| cross-check | `kron` |
| support | `examples` (on-demand real data), `gaiadb` (the opt-in offline Gaia catalogue), `astap` (finding/fetching the solver) |

## Three FITS layouts, two header dialects

All normalise to the same in-memory frame. See `docs/data-format.md`.

| | `"cube"` — native stack | `"mef"` — CrowdSky | `"bayer"` — raw sub |
|---|---|---|---|
| structure | 3-D primary HDU, FITS axes `(nx, ny, 3)` → numpy `(3, ny, nx)` | empty primary + `RED`/`GREEN`/`BLUE` ImageHDUs (+ `FOOTPRINT`, `STAR-TAB`) | single 2-D mosaic, `BAYERPAT` (`GRBG`); demosaiced on load, native samples kept on `frame.bayer` |
| sub count | `STACKCNT` | `NIMAGES` | — |
| `EXPTIME` | **per sub-exposure** (10 s) | **total on-sky** (e.g. 410 s) | per sub |
| total on-sky | `TOTALEXP` | `EXPTIME` | — |
| exposure span | `DATE-OBS` only | `OB-START` / `OB-END` (use these) | `DATE-OBS` |

That `EXPTIME` collision is the trap: reading it naively gives a per-sub time of
410 s and a mid-exposure timestamp ~5 minutes wrong. `frames.frame_metadata`
resolves the dialect from *which keywords are present*, not from the layout; never read
these keywords directly. `stacking` writes its output in the **native** dialect so a local
stack reads back exactly like an on-board one (`layout="stacked"`).

## The batch runner

`pipeline` has three stages (`solve_all` → `build_frame_table` → `build_measurements`),
run in that order, over one shared `_run`. Its contract, which new stages must keep:

- `work_fn(key, project)` returns `(key, status, payload)` and **never raises** — a raise
  inside a worker loses the frame identity that makes the message actionable;
- statuses are `ok / cached / no_wcs / load_error / failed`; one bad frame never sinks
  the batch;
- **idempotent and checkpointed** — killing a stage and restarting it is a normal
  operation, not a recovery procedure;
- the catalogue is loaded **once per worker** via the pool initializer;
- threads instead of processes when `workers == 1`, and capped at 4 for `solver="nova"`;
- retries on *network* stages only. A failed SEP extraction is not transient.

`FrameSource` (`keys()` / `path(key)`) is the seam for a future remote archive.
Implementations must be picklable — the runner ships them to worker processes.

## Tests

- `tests/conftest.py` injects Gaussian PSFs of **known** flux and builds the reference
  catalogue by inverting the exact calibration relation the pipeline fits, so tests assert
  *recovery* of `ZP_TRUE` / `K_TRUE` / aperture / period, not merely "it ran".
  If you change the calibration relation, `expected_zeropoint()` is where the pinning lives.
- `tests/test_real_data.py` runs against real MW Cam cutouts fetched by
  `examples.download()` (~18 MB, cached; `SEESTAR_PHOTOMETRY_DATA` overrides the location).
  The module skips itself when offline — a green run on a disconnected machine is not proof.
- Every documented example runs on that same `examples` data, on purpose: an example that
  can't be executed rots.

## Style

- numpydoc docstrings, rendered by Sphinx napoleon. Comments explain **why**, and the
  measured numbers behind a choice go in the docstring — that density is the house style,
  not clutter.
- Source stays ASCII: `--` in docstrings, not an em dash. Markdown files may use one.
- Lines wrap under ~95 characters.
- Module-level constants are documented with `#:` so autodoc picks them up.

## Environment

`ASTROMETRY_KEY` (astrometry.net) is read from the environment and never stored in
the repo — only needed for `solver="nova"`; the default ASTAP solver is offline and
expects the binary at `astrometry.ASTAP_EXE`.

## Where the "why" lives

`docs/` are decision records, not API reference. Read the relevant one before changing a
numerical choice: `photometry-design.md` (apertures — start here), `data-format.md`,
`astrometry-and-gaia.md`, `frame-table.md`, `light-curves.md`, `diagnostics.md`,
`architecture.md`. Update `CHANGELOG.md` for anything user-visible.
