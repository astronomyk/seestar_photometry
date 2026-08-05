# How deep does it go?

Two questions, and they need different answers: *how deep is this frame*, and *how deep would
N hours be*.

## This frame

```python
from seestar_photometry import calibration, examples, photometry, quality

frame = examples.stack()
ext = photometry.extract_sources(frame)
ext.match_gaia(examples.gaia(), wcs=examples.wcs())
cal = calibration.fit_zeropoint(ext.sources, band="G")

row = quality.frame_quality(ext, cal)
print(f"{row['total_exptime']:.0f} s  ->  V = {row['v_lim_5sigma']:.2f} at 5 sigma")
# 760 s  ->  V = 17.15 at 5 sigma
```

The limit comes from the zero point and the aperture noise, not from counting faint
detections:

$$V_\mathrm{lim} = \mathrm{ZP} - 2.5\log_{10}(n\,\sigma_\mathrm{aper})$$

where `sigma_aper` is the aperture flux error — roughly √(aperture area) × sky RMS. That
means it needs no detections near the limit at all; it extrapolates the *bright*-star
calibration through the noise model.

```{admonition} This number is optimistic by about 0.7 mag
:class: warning
`sigma_aper` assumes independent pixels. Demosaic interpolation and sub-pixel resampling
**correlate** neighbours, and the measured aperture noise on S50 on-board stacks is **2.0×**
the independent-pixel prediction.

Relative use is fine — the offset is near-constant within one processing chain, so
exposure-time fits and per-unit rankings stand. Absolute depth claims are not, and the value
must never be compared across processing chains. An earlier validation against the faint edge
of the V–SNR locus looked reassuring but was circular: both sides used the same per-pixel
error model, so it could not detect an error *in* that model.
```

## Longer exposures

Background-limited photometry deepens by 1.25 mag per decade of exposure. `depth` applies
that to put frames of different integration on a common footing:

```python
from astropy.table import Table
from seestar_photometry import depth

frames_table = Table(rows)        # one row per frame, from pipeline.build_frame_table

print(depth.detection_limit(frames_table, t_ref=900.0))
# model   exptime  n_frames  v_lim_median  v_lim_std
# S50       900.0         3         17.25       ...

print(depth.limit_vs_exptime(frames_table, exptimes=(900, 3600, 14400)))
# 900 s -> 17.25    3600 s -> 18.00    14400 s -> 18.75
```

So a 15-minute S50 stack reaches V ≈ 17.25 and four hours would reach ≈ 18.75 — *if* √t holds
that far.

## Does √t actually hold?

Test it rather than assuming. The two bundled c17 stacks are the same field and unit at two
integrations:

| Exposure | 5σ limit |
|---|---|
| 760 s | 17.15 |
| 1460 s | 17.69 |

Measured gain **+0.53 mag**; the √t law predicts **+0.35 mag**. The real frame did better
than the ideal scaling, because the deeper stack also happened to have slightly better seeing
and a darker sky — which is exactly why exposure alone is not the whole story.

## Conditions, not just time

Depth is set by integration *and* conditions. `fit_depth_model` separates them:

$$V_\mathrm{lim} = a + b\log_{10} t + c\,\mathrm{SQM} + d\,\mathrm{FWHM}$$

```python
model = depth.fit_depth_model(frames_table)
print(model["coeffs"], model["rms"], model["n_frames"])
```

Fit it **per unit** — each Seestar has its own response. Three things to know before using it:

- **It is a noise model, not a zero-point correction.** It acts on the faint-end
  sensitivity, never on the flux scale, so it does *not* change a detected star's magnitude
  or its error bar.
- **It mostly re-derives photon shot noise.** For a sky-limited stack with an aperture ∝ FWHM,
  the expected coefficients are ≈ +0.5 per mag/arcsec² of sky and ≈ −0.3 per pixel of FWHM.
  A real fit landed at +0.46 and −0.20. The instrument-specific signal is in the
  *deviations*, and in the differences *between* units.
- **Beware a narrow range.** A unit whose frames span little variation in sky or seeing has
  badly constrained `c` and `d`, and they will happily take implausible values.

What it is genuinely good for: predicting the depth of a planned observation, setting 5σ
upper limits for non-detections, weighting frames in a co-add, and flagging a frame that
falls far below its condition-predicted depth (cloud, trailing, focus).

One real lever it suggests: since FWHM drives noise through aperture *area*, a tighter
aperture plus an aperture correction can claw back some depth in poor seeing.

## Looking at it

```python
from seestar_photometry import plots

plots.depth_vs_exptime(frames_table)             # with the sqrt(t) line drawn
plots.depth_vs_driver(frames_table, "sky_sb")    # exposure divided out
```

Flattening at the long-exposure end in the first panel means a systematic floor has taken
over and the extrapolation no longer applies. Validate a multi-hour claim by actually
co-adding frames — see [stacking](stacking.md) — rather than by extending the line.

Full detail in [frame-table](../frame-table.md).
