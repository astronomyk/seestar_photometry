# Migrating from `mwcam` / `mw-cam-py`

This package supersedes the `mwcam` engine in the `photometry_adventures_2` repo, which
stays frozen there for paper reproducibility. Nothing in that repo was changed.

## Why

The algorithms were settled and documented; the *packaging* blocked reuse:

1. **Three live copies of the engine** — `src/mwcam/`, and two vendored `_phot/` trees. The
   supernova copy was byte-identical, but `mwcam_catalogue/_phot/calibration.py` had
   **forward-evolved** (300 vs 225 lines) and the "real" package never received those
   improvements. Fixes landed in whichever copy happened to be open.
2. **Every dataset-specific constant was a module global** in a driver script:
   `DEFAULT_DATA_DIR`, `GAIA_DIR`, `EXTRA_DIRS`, `ENCLOSED`, `FIELD`, `FIT_MAG_RANGE`,
   `DATA_ROOT`, `PAPER_DATA`, `ASTAP_EXE`, `SN_RA_DEG`. Retargeting meant editing code.
3. **The batch driver was copy-pasted four times** (`solve_all_wcs.py`,
   `build_frame_quality.py`, `build_lightcurve_table.py`, `mwcam_catalogue/catalogue.py`).
4. **Two names for one table**, `frame_quality.ecsv` and `frame_catalogue.ecsv`, with
   diverging column sets.
5. **No tests at all** — including for the lazy-matplotlib rule the docs claimed was verified.

## Module map

| `mwcam` | `seestar_photometry` | Notes |
|---|---|---|
| `io` | `frames` | Extended: both FITS layouts, both header dialects, `FrameSource`/`LocalTree` |
| `photometry` | `photometry` | Same algorithms; sizing factored into `_per_band_setup` so detection and forced photometry cannot drift apart; `Extraction` now keeps `cogs` for the figures |
| `astrometry` | `astrometry` | ASTAP and astrometry.net unified behind `solve()`; `lift()` added |
| `catalogs` | `catalogs` | Retry-with-cache-on-success-only folded in from `catalogue.py` |
| `calibration` | `calibration` + `quality` | The **evolved** fork adopted; `frame_quality` split out |
| `variability` | `lightcurves` | Plus `build_stars`, `comparison_curves`, `periodogram`, `phase_fold` |
| `detection` | `depth` | Plus `fit_depth_model` |
| `diagnostics` | `kron` | Renamed — "diagnostics" now means the figure layer |
| — | `project`, `pipeline` | New: config and the unified runner |
| — | `plots`, `report`, `_style` | New: the figure layer |
| — | `contamination` | Extracted from the supernova module |

## Behaviour changes

**The evolved calibration is now canonical.** From `mwcam_catalogue/_phot/calibration.py`:
`v_jkc_mag → v_mag` fallback, the `b_jkc−r_jkc → bp_rp → none` colour cascade,
`saturation_mag`/`v_sat`, and `v_lim_100sigma`.

**Two aperture fractions, explicitly.** `mwcam` documented `enclosed=0.90` but
`build_lightcurve_table.py` silently used `0.95`. Both are right for their purpose, so
`Project` carries `enclosed_characterise=0.90` and `enclosed_lightcurve=0.95`, with the
reason recorded. See [photometry-design.md](photometry-design.md).

**ASTAP is the default solver**, not astrometry.net — local, offline, no key, no rate limit.

**Timing is fixed for CrowdSky frames.** `mwcam.io.frame_metadata` read `EXPTIME` at face
value, which is the *total* on a CrowdSky frame, not per-sub. `total_exptime` came back `nan`,
so `frame_times` used a half-exposure of 0 and the mid-exposure timestamp ran **~5 minutes
early** — enough to smear a transit. `frames._exposure` resolves the dialect, and
`frame_times` prefers the true `OB-START`/`OB-END` span when present.

**`frame_quality.ecsv` + `frame_catalogue.ecsv` → `frames.ecsv`**, one superset schema. Extra
project columns go through `Project.provenance`.

**Cross-match column names**: `gaia_ra`/`gaia_dec` → `cat_ra`/`cat_dec`, since the reference
catalogue need not be Gaia.

**Frame curation is not carried over.** `io.find_frames(max_flat=18)` hard-coded an
MW-Cam-specific rule (dithers `flat_19`–`flat_22` were tree-obstructed). `LocalTree` takes an
optional `curate(path) -> bool` predicate instead; no rule is applied by default. Add the
dither rule as a predicate if and when you need it.

## Porting a driver script

`build_frame_quality.py` (115 lines of globals + a ProcessPool) becomes:

```python
from seestar_photometry import Project, Target, LocalTree, pipeline

proj = Project(
    target   = Target("MW Cam", ra=186.6821, dec=81.474),
    source   = LocalTree(roots=[r"D:\seestar_paper_data_sets_2\MW Cam s50\stacks"]),
    work_dir = r"D:\work\mwcam",
)
pipeline.solve_all(proj)
pipeline.build_frame_table(proj, diagnostics=3)
```

The `EXTRA_DIRS` pattern becomes extra entries in `LocalTree(roots=[...])`.

## Reproducing the old numbers

To compare against a frozen `frame_catalogue.ecsv`, point a `Project` at the same stacks with
`enclosed_characterise=0.90` and reuse the existing `.wcs` sidecars — they live beside the
frames, so astrometry is held fixed and only the photometry is compared. `examples/` includes
the MW Cam driver used for exactly this.
