# seestar-photometry — conventions

Read this before changing code. The rationale for every numerical choice lives in
`docs/`; this file is the shape of the package.

## Names

Distribution `seestar-photometry`, import package `seestar_photometry`. Python
import names can't contain hyphens.

## Hard conventions

- **`BANDS = ("R", "G", "B")` is the canonical axis-0 order.** Every frame is a
  `(3, ny, nx)` float32 array in that order. Every per-band array (`Extraction.rms`,
  `.aperture`, `.fwhm`) is shape `(3,)` in that order. **Green is index 1 and is the
  science band** (≈ Johnson V through the Seestar's IRCUT filter).

- **Image shape comes from the data, never the header.** `frame.g.shape` gives
  `(ny, nx)`. The CrowdSky multi-extension layout has an *empty* primary HDU with no
  `NAXIS1`/`NAXIS2`, so any code reading the shape from a header is broken on half
  the datasets.

- **`SeestarFrame` and `Project` are thin data containers, not god objects.** All
  analysis is module-level *functions* taking the object as their first argument
  (`photometry.extract_sources(frame, ...)`). Methods that briefly lived on the
  frame class were moved out on purpose.

- **Optional dependencies are imported lazily**, inside the function that needs
  them — never at module top. `matplotlib` in particular: `import
  seestar_photometry` must work in a core install. `tests/test_imports.py` asserts
  `matplotlib` is absent from `sys.modules` after import, so this can't rot.

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
| data | `frames` |
| measure | `photometry`, `astrometry`, `catalogs` |
| calibrate | `calibration`, `quality` |
| science | `lightcurves`, `contamination`, `depth` |
| orchestrate | `project`, `pipeline` |
| inspect | `plots`, `report` |
| cross-check | `kron` |

## Two FITS layouts, two header dialects

Both normalise to the same in-memory frame. See `docs/data-format.md`.

| | `"cube"` — native Seestar | `"mef"` — CrowdSky |
|---|---|---|
| structure | 3-D primary HDU, FITS axes `(nx, ny, 3)` → numpy `(3, ny, nx)` | empty primary + `RED`/`GREEN`/`BLUE` ImageHDUs (+ `FOOTPRINT`, `STAR-TAB`) |
| sub count | `STACKCNT` | `NIMAGES` |
| `EXPTIME` | **per sub-exposure** (10 s) | **total on-sky** (e.g. 410 s) |
| total on-sky | `TOTALEXP` | `EXPTIME` |
| exposure span | `DATE-OBS` only | `OB-START` / `OB-END` (use these) |

That `EXPTIME` collision is the trap: reading it naively gives a per-sub time of
410 s and a mid-exposure timestamp ~5 minutes wrong. `frames.frame_metadata`
resolves the dialect; never read these keywords directly.

## Environment

`uv` only. `uv sync --extra dev`, then `uv run pytest` / `uv run python ...`.
Never `uv run --active`. Point PyCharm at `.venv`; don't let it make its own.

`ASTROMETRY_KEY` (astrometry.net) is read from the environment and never stored in
the repo — only needed for `solver="nova"`; the default ASTAP solver is offline.
