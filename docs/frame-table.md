# The frame table

`quality.frame_quality(extraction, cal)` returns one row of per-frame metrics;
`pipeline.build_frame_table` accumulates them into **`work_dir/frames.ecsv`**, one row per
frame. It runs process-parallel — WCS and the catalogue are cached, so it is pure CPU — and
never touches the original FITS.

This replaces two diverging tables in the predecessor (`frame_quality.ecsv` and
`frame_catalogue.ecsv`) with one schema. Project-specific bookkeeping goes through the
`Project.provenance` hook rather than a second builder.

## Columns

| Group | Columns |
|---|---|
| identity | `path`, `frame`, `layout`, `model`, `unit`, `telescope`, `object`, `filter` |
| timing / geometry | `date_obs`, `obs_start`, `obs_end`, `eqmode`, `n_exp`, `exptime`, `total_exptime`, `airmass`, `site_lat`, `site_lon`, `ccd_temp`, `pixscale` |
| counts | `n_sources`, `n_green`, `n_matched`, `n_cal` |
| calibration | `zeropoint(_err)`, `colour_term(_err)`, `colour0`, `rms`, `chi2_red`, `fit_mag_lo/hi` |
| depth / range | `v_lim_5sigma`, `v_lim_100sigma`, `v_sat`, `sigma_aper` |
| PSF | `fwhm_R/G/B`, `aperture_R/G/B` |
| sky | `sky_R/G/B`, `sky_pedestal`, `bg_poly`, `bg_resid_std`, `sky_sb`, `bortle` |
| astrometry | `median_arcsec`, `p90_arcsec`, `matched_frac`, `n` |
| provenance | `solver`, `enclosed`, `onboard_*` (CrowdSky only), plus anything from the hook |

`eqmode` is 0 for Alt-Az and 1 for equatorial, and is absent from CrowdSky headers.
**`airmass` is computed**, not read — Seestar headers carry no `AIRMASS`.

## Reading the metrics

**`rms`** — green-band calibration scatter over the fit window. The master quality number.

| `rms` | verdict |
|---|---|
| 0.02–0.04 | photometric |
| < 0.06 | usable (the conventional cut, `quality.photometric`) |
| 0.1–0.3 | marginal |
| > 0.5 | junk — cloud, trailing, bad focus |

**`chi2_red`** — scatter relative to the photon-noise prediction (`1.0857/SNR`). This is a
*relative* indicator, **not** an absolute "≈ 1" target: even the best frames sit at 100–200,
because the ~0.03 mag floor is **systematic** (flat-field, intra-frame PSF), not photon
noise. Junk frames reach 1e4–1e5. Use it to find frames much worse than their peers.

**`n_cal`** — a frame whose calibration-star count collapses should be distrusted whatever
its `rms` says: with few stars the sigma clipping has nothing to work with, so the scatter
can look deceptively small.

**`sky_sb` / `bortle`** — V-equivalent sky surface brightness,
`ZP − 2.5·log₁₀((pedestal − BIAS) / pixscale²)`. It carries a sensor- and band-dependent
offset (V calibrated through the broad IRCUT band), so compare *within* a unit rather than
blindly across models — on one site the S50 read ~19.8 and the S30pro ~18.8 mag/arcsec².

**`bg_poly`** — the six `photometry.fit_background` coefficients (`1, x, y, x², xy, y²`);
`bg_resid_std` is the residual after that smooth fit (small = smooth sky, large =
bright-star halos).

**`v_sat`** — the empirical bright-end limit: the faintest matched star whose peak pixel is
at or above 95% of the 16-bit ceiling. `nan` means nothing in the field saturates.

## `v_lim_5sigma` — and why it is optimistic

`calibration.limiting_mag` computes `V_lim = ZP − 2.5·log₁₀(5·σ_aper)`, where `σ_aper` is
the constant aperture flux error (median green `fluxerr` ≈ √(aperture area) × sky RMS). It
needs no detections near the limit — it extrapolates the bright-star zero point through the
noise model.

> **It is optimistic in absolute terms by ~0.7 mag.** `σ_aper` assumes independent pixels,
> but demosaic interpolation and sub-pixel resampling **correlate neighbours**: the measured
> aperture noise is **2.0× the independent-pixel prediction** on S50 on-board stacks.
>
> Relative use is fine — the offset is near-constant within one processing chain, so
> exposure-time fits and per-unit rankings stand. Absolute depth claims are not, and the
> column must **never** be compared across processing chains.
>
> An earlier validation against the faint edge of the V–SNR locus (~16.8 analytic vs ~16.6
> measured) looked reassuring but was largely circular: both sides use the same per-pixel
> error model, so it could not detect an error *in* that model.

The S30pro runs ~1–1.8 mag shallower than the S50s, as expected from its smaller aperture.

## Depth vs exposure (`depth`)

`v_lim_5sigma` is per frame, at that frame's own `total_exptime` (a non-round value — a
stack is `n_exp` × 10–20 s). `depth` rescales it via the background-limited √t law:

    V_lim(t) = V_lim(t₀) + 1.25·log₁₀(t/t₀)

`depth.detection_limit(frames)` gives the 15-minute 5σ limit per group (S50 ≈ 17.2,
S30pro ≈ 16.7); `depth.limit_vs_exptime(frames)` traces 1 min – 16 h. The long end assumes
√t still holds, i.e. no systematic floor — validate by actually co-adding frames before
quoting it.

## Condition-corrected depth — what it is, and what it is *not*

The scatter in `v_lim` vs exposure shrinks markedly once per-frame conditions are modelled
(`depth.fit_depth_model`):

    v_lim = a + b·log₁₀(t) + c·SQM + d·FWHM

Fit **per unit** — each Seestar gets its own `c`, `d`.

- **It is a noise/depth model, not a zero-point correction.** It acts on the aperture-noise
  term (faint-end *sensitivity*), never on the flux scale. The flux→magnitude calibration is
  the per-frame zero point, fit from bright stars and already independent of sky conditions.
  So this correction **does not change a detected star's magnitude or its error bar**.
- **It is largely just the background-limited noise law.** For a sky-limited stack with an
  aperture ∝ FWHM, noise ∝ FWHM and sky shot noise ∝ 10^(−0.2·SQM), so the expected
  coefficients are **≈ +0.5 per mag/arcsec² of SQM** and **≈ −0.3 per pixel of FWHM**. The
  c17 fit (c ≈ +0.46, d ≈ −0.20) sits right on the SQM prediction and in the same ballpark
  on FWHM — the correction mostly *re-derives photon-shot-noise physics*. The genuinely
  instrument-specific signal is in the **deviations** and in the **between-unit differences**
  in `c`, `d`.
- **Beware narrow-range units.** A unit whose frames span only a small range of sky or
  seeing has poorly constrained `c` and `d`, which will happily take implausible values.

What it is actually good for — as a predictor *across* frames, not a correction *within* one:

- exposure-time planning: invert it to predict the depth of a planned observation;
- principled 5σ upper limits for non-detections;
- frame weighting in co-adds (the one place it improves combined photometry);
- quality control: a frame far below its condition-predicted `v_lim` is suspect;
- characterisation: how conditions cap the science yield.

One real photometric lever it *suggests*: because FWHM drives noise through aperture *area*,
a tighter aperture plus an aperture correction can claw back some depth in poor seeing.

## Cross-unit finding

The S30pro came out markedly more photometrically stable (median `rms` ~0.05) than the three
S50 units (~0.12–0.14) — a genuine cross-unit difference, not a processing artefact.

## Re-running as a dataset grows

Everything is idempotent, so after adding frames just re-run the stages in order:

```python
pipeline.solve_all(proj)          # solves only the new frames
pipeline.build_frame_table(proj)  # appends only the new rows
```

A degenerate, non-celestial solution is rejected rather than cached, and astrometry.net
throttles above ~4 concurrent solves — see [astrometry-and-gaia.md](astrometry-and-gaia.md).
