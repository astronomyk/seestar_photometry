# seestar-photometry

Time-domain photometry from ZWO Seestar smart-telescope stacks, calibrated onto Gaia
DR3 synthetic Johnson V.

Built from two working pipelines — a δ-Scuti light curve of MW Cam (P = 0.1294 d, a
~23 mmag scatter floor) and a Type II supernova in NGC 3310 — generalised so the same
recipe applies to a new target with a short driver script and no code edits.

The recipe:

1. load a stack (either FITS layout — native Seestar or CrowdSky);
2. detect sources with SEP and size a circular aperture **per frame and per band** from
   the curve of growth;
3. solve a per-frame WCS and cache it beside the FITS — the on-board WCS is off by
   ~1 arcmin and unusable for photometry;
4. cross-match a once-cached Gaia mosaic carrying synthetic JKC V;
5. fit `V = m_inst + ZP + k·(B−R)` on the green plane;
6. forced-aperture photometry at fixed sky positions, so no series ever goes ragged;
7. per-frame ensemble zero point from N comparison stars, **each referenced to its own
   catalogue magnitude**.

## Install

```bash
uv sync --extra dev
```

Python ≥ 3.11. Core deps are `astropy`, `astroquery`, `numpy`, `scipy`, `sep`;
`matplotlib` is the `plot` extra (the measurement path runs headless).

For plate solving, either install [ASTAP](https://www.hnsky.org/astap.htm) (local,
offline, the default) or set `ASTROMETRY_KEY` for astrometry.net.

## Quickstart

```python
from seestar_photometry import Project, Target, LocalTree, pipeline, lightcurves, report

proj = Project(
    target   = Target("MW Cam", ra=186.6821, dec=81.474),
    source   = LocalTree(roots=[r"D:\data\MW Cam s50\stacks"]),
    work_dir = r"D:\work\mwcam",
)

pipeline.solve_all(proj)                                  # .wcs sidecars, idempotent
frames = pipeline.build_frame_table(proj, diagnostics=3)   # frames.ecsv + figures
stars, meas = pipeline.build_measurements(proj)            # stars.ecsv + measurements.ecsv

comps = lightcurves.select_comparisons(stars, dmag=1.0, colour_tol=0.3,
                                       max_sep_arcmin=15)
lc = lightcurves.differential_lightcurve(
    meas, lightcurves.target_id_of(stars), comps, band="G"
)
report.lightcurve_report(lc, stars, meas, comps, proj.diagnostics_dir)

print(f"{len(lc)} epochs, scatter {lc.meta['scatter'] * 1000:.0f} mmag")
```

All three stages are resumable — interrupt and re-run, they pick up where they stopped.
Run them in order: stages 2 and 3 read the cached WCS and never solve.

See `examples/` for complete drivers, including a supernova with host-galaxy
subtraction and a template for a new dataset.

## Outputs

Everything derived lands in `work_dir` (never in the data tree, never in FITS headers).
The one exception is the per-frame `.wcs` sidecar, which lives beside its frame because
it is expensive to recompute and useful to every project touching that frame.

| File | What it is |
|---|---|
| `frames.ecsv` | one row per frame: zero point, colour term, scatter, PSF, sky, depth limits |
| `stars.ecsv` | one row per catalogue source measured, with separation from the target |
| `measurements.ecsv` | the long table: one row per (source, frame, band) |
| `diagnostics/*.png` | the figure sets below |

## Diagnostics

Pass `diagnostics=` to a pipeline stage, or call `report.*` directly. Three sets:

- **Per frame** — curve of growth with the chosen aperture marked, per-band FWHM and
  its chromatic ratios, the zero-point relation and colour-term fit, residuals against
  magnitude / SNR / radius, a residual map over the frame, the background triptych,
  detections with apertures drawn on, and the cross-match separation histogram.
- **Per dataset** — zero point and conditions over time, `rms` and `chi2_red`
  distributions against the photometric-grade cut, depth against exposure with the
  √t law overlaid and against its condition drivers, calibration coverage, and a
  one-page contact sheet.
- **Per light curve** — finder chart, **each comparison's own differential curve**
  (the single most informative check — a variable or blended comparison shows up
  immediately, sorted worst-first), the ensemble zero point with its star-to-star
  spread, achieved scatter against comparison brightness, periodogram, and phase fold.

## Documentation

`docs/` holds the decision records — *why* each numerical choice is what it is, with
the measurements behind it. Start with `docs/photometry-design.md`.
`CLAUDE.md` holds the conventions to follow when changing code.

## Tests

```bash
uv run pytest
```

Fully offline: synthetic frames with injected Gaussian PSFs of known flux, in both FITS
layouts, so the tests assert *recovery* of a known zero point, colour term, aperture and
period rather than merely that the code runs.
