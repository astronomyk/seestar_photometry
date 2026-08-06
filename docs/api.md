# API reference

Modules are layered: each imports only from the ones above it. The docstrings carry the
detail — this page is the map.

| Layer | Module | Role |
|---|---|---|
| data | `frames` | Load any layout into one `(3, ny, nx)` cube; resolve header dialects; discover frames |
| data | `debayer` | Bayer demosaic for raw subs |
| measure | `photometry` | **The standard path**: detection, per-band aperture, fixed and forced flux |
| measure | `astrometry` | Per-frame WCS against the catalogue, ASTAP or astrometry.net, cached as a sidecar |
| measure | `catalogs` | Cached Gaia mosaic, footprint subsetting, cross-match |
| measure | `stacking` | Register and co-add raw subs |
| support | `gaiadb` | The offline Gaia catalogue: an opt-in local copy of the TAP query |
| support | `astap` | Locating the ASTAP solver, and fetching one when there is none |
| calibrate | `calibration` | Zero point + colour term; limiting and saturation magnitudes |
| calibrate | `quality` | The frame-table row |
| science | `lightcurves` | Timing, comparison selection, ensemble differential, periods |
| science | `contamination` | Host-galaxy subtraction |
| science | `depth` | √t rescaling, detection limits, the condition model |
| orchestrate | `project` | `Target`, `Project` — the whole configuration surface |
| orchestrate | `pipeline` | The three batch stages over one resumable runner |
| inspect | `plots`, `report` | Diagnostic figures and saved figure sets |
| cross-check | `kron` | Kron/AUTO photometry — never the science path |
| — | `examples` | The downloadable real example data |

Three conventions hold everywhere:

- `BANDS = ("R", "G", "B")` is the canonical axis-0 order, and **green is index 1 and is the
  science band**. Every per-band array is shape `(3,)` in that order.
- **Image shape comes from the data, never a header.** The CrowdSky layout has an empty
  primary HDU with no `NAXIS1`/`NAXIS2`.
- `SeestarFrame` and `Project` are thin data containers. Analysis lives in module-level
  functions taking them as the first argument.

## Data

```{eval-rst}
.. automodule:: seestar_photometry.frames
.. automodule:: seestar_photometry.debayer
.. automodule:: seestar_photometry.examples
.. automodule:: seestar_photometry.gaiadb
.. automodule:: seestar_photometry.astap
```

## Measurement

```{eval-rst}
.. automodule:: seestar_photometry.photometry
.. automodule:: seestar_photometry.astrometry
.. automodule:: seestar_photometry.catalogs
.. automodule:: seestar_photometry.stacking
```

## Calibration

```{eval-rst}
.. automodule:: seestar_photometry.calibration
.. automodule:: seestar_photometry.quality
```

## Science

```{eval-rst}
.. automodule:: seestar_photometry.lightcurves
.. automodule:: seestar_photometry.contamination
.. automodule:: seestar_photometry.depth
```

## Orchestration

```{eval-rst}
.. automodule:: seestar_photometry.project
.. automodule:: seestar_photometry.pipeline
```

## Diagnostics

```{eval-rst}
.. automodule:: seestar_photometry.plots
.. automodule:: seestar_photometry.report
```

## Cross-checks

```{eval-rst}
.. automodule:: seestar_photometry.kron
```
