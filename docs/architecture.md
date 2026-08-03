# Package architecture

## Names and environment

**Distribution `seestar-photometry`, import package `seestar_photometry`.** `src/` layout,
hatchling backend, installed editable into a **uv** environment.

Core deps: `astropy`, `astroquery`, `numpy`, `scipy`, `sep`. `matplotlib` is the `plot`
extra so the measurement path installs and runs headless. Deliberately *not* dependencies:
`photutils` (unused — all photometry is `sep`), `seestarpy` (declared but never imported by
the predecessor), `crowdsky-client` (a future optional extra).

`uv run pytest`, `uv run python ...`. Never `uv run --active`.

## Module layers

Each module imports only from layers above it.

| Layer | Module | Role |
|---|---|---|
| data | `frames` | Load either FITS layout into one `(3, ny, nx)` frame; resolve header dialects; discover frames (`FrameSource`, `LocalTree`) |
| measure | `photometry` | **The standard path**: detection, per-band aperture sizing, fixed-aperture and forced flux, curve of growth, background fit |
| measure | `astrometry` | Per-frame WCS via ASTAP or astrometry.net, cached as a `.wcs` sidecar |
| measure | `catalogs` | Cached oversized Gaia mosaic, footprint subsetting, cross-match |
| calibrate | `calibration` | Zero point + colour term against synthetic V; limiting and saturation magnitudes; sky brightness |
| calibrate | `quality` | The `frames.ecsv` row builder and schema |
| science | `lightcurves` | Timing, comparison selection, ensemble differential, period tools |
| science | `contamination` | Extended-emission (host galaxy) subtraction |
| science | `depth` | √t rescaling, detection limits, the condition-corrected depth model |
| orchestrate | `project` | `Target`, `Project` — the whole configuration surface |
| orchestrate | `pipeline` | The three batch stages over one resumable runner |
| inspect | `plots`, `report` | Diagnostic figures, and saved figure sets |
| cross-check | `kron` | Kron/AUTO photometry — never the science path |

## Conventions, and why

**`BANDS = ("R", "G", "B")` is the canonical axis-0 order.** Every per-band array
(`Extraction.rms`, `.aperture`, `.fwhm`) is shape `(3,)` in that order. **Green is index 1
and is the science band** (≈ Johnson V through IRCUT).

**Image shape comes from the data, never a header.** The CrowdSky layout has an *empty*
primary HDU with no `NAXIS1`/`NAXIS2`, so any code reading shape from a header is broken on
half the datasets. `frame.shape` is the accessor.

**`SeestarFrame` and `Project` are thin data containers, not god objects.** All analysis is
module-level *functions* taking the object as their first argument
(`photometry.extract_sources(frame, ...)`). Methods that briefly lived on the frame class in
the predecessor were moved out for exactly this reason.

**Configuration is Python, not files.** No TOML, no console scripts. Everything the
predecessor kept in module-level globals — data roots, target coordinates, aperture
fractions, fit windows, solver paths — lives on `Project`. Retargeting means constructing a
different `Project`, never editing code.

**Optional dependencies are imported lazily**, inside the function that needs them.
`matplotlib` in particular: `import seestar_photometry` must work in a core install. The
predecessor *documented* this rule but never tested it;
`tests/test_imports.py` now asserts `matplotlib` (and `astroquery`) are absent from
`sys.modules` after import, in a fresh interpreter so the test can't fool itself.

**The standard path and the cross-checks are separate modules.** `photometry` is what the
science uses; `kron` exists for star/galaxy separation and total-flux comparison. Keeping
Kron out of `photometry` makes the science path unambiguous.

**Derived products never go into the data tree or the FITS headers.** Tables and figures go
to `Project.work_dir`. The single exception is the per-frame `.wcs` sidecar, which lives
beside its frame because it is expensive to recompute, is a property of the frame rather than
of any analysis, and is reused by every project touching that frame.

**No defensive error handling in the measurement layer.** A broken frame raises. The batch
runner catches per-frame exceptions, records a status, logs it, and continues. Keeping those
two concerns apart is why the science code reads cleanly.

## The batch runner

`pipeline._run` replaces four copy-pasted drivers in the predecessor. One contract:

- `ProcessPoolExecutor` with an initializer that loads the catalogue **once per worker** —
  the difference between seconds and minutes over a large dataset;
- **threads instead of processes when `workers == 1`**, since a single worker gains nothing
  from a separate process and pays the full spawn cost (on Windows, re-importing astropy and
  scipy before any work starts);
- threads capped at 4 for `solver="nova"`, because astrometry.net throttles above that;
- per-frame `try/except` → `(key, status, payload)` with statuses
  `ok / cached / no_wcs / load_error / failed`; one bad frame never sinks the batch;
- **idempotent** — a frame already in the output table, or already carrying a sidecar, is
  skipped unless `force=True`;
- **checkpointed** every N frames, so a killed run loses at most that many;
- retries on *network* stages only, never re-running SEP: a failed extraction is not
  transient, and retrying it multiplies the cost of a frame that was never going to work.

Killing a stage and restarting it is a normal operation, not a recovery procedure — which
matters when a stage takes hours over thousands of frames.

## The frame-source seam

`FrameSource` is a two-method protocol: `keys()` and `path(key)`. That is all the pipeline
needs, and it exists so a remote archive — one that discovers frames through an API and
downloads them on demand — can be added later without touching pipeline code. Only
`LocalTree` is implemented today.

Implementations must be picklable, since the runner ships them to worker processes.
