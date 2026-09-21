# Passbands for the S30 Pro and S50 Pro

The S50's photometric system is characterised in
[the S50 photometric system](s50-photometric-system.md). This page does the same two
newer models, by the same method, and records the one place that method had to change —
along with the control that says the change is real rather than convenient.

Both units observed the same field, 20 Cep, as the S50 did: 16 stacks from two S30 Pro
units and 24 from three S50 Pro units, all from `D:\seestar_saturation_data_set`. Built by
`tools/digitise_sensor_qe.py` and `tools/fit_unit_passbands.py`.

## What is fitted, and against what

A passband here is the sensor's quantum efficiency times a UV/IR-cut filter model. The QE
comes from the manufacturer's plot, digitised by colour-masking; the filter is never
published for any of these units, so it is **fitted from the photometry**. Calibrate each
Bayer plane against a synthetic magnitude computed at a trial filter, measure the slope of
the zero-point residual against colour, and take the value where that slope crosses zero.

The colour axis is Gaia BP−RP, and on these units it is **synthesised from the same XP
spectra** as the reference magnitudes, because `bp_rp` lives only in `gaia_source` and is
masked in the offline all-sky catalogue. Checked against real Gaia colours for 12,948
stars in this field, the synthetic colour agrees to **0.0102 mag** scatter about a
constant −0.351 offset, which is just AB−Vega. A constant cannot move a zero crossing, and
neither can a scale error, so the fitted edge is insensitive to both.

## The single-edge model breaks on both new sensors

For the S50, one IR-cut edge satisfied the green and red planes simultaneously: they
crossed zero 1 nm apart. On both new units they do not.

| solo crossing | G | R | gap |
|---|---|---|---|
| S50 | 693.3 nm | 692.8 nm | **0.5 nm** |
| S30 Pro | 703.2 nm | 678.4 nm | **24.8 nm** |
| S50 Pro | 708.5 nm | 684.9 nm | **23.5 nm** |

There is one filter, so both planes must be fitted by one edge — and at the compromise
value (their mean) each plane is left with a residual colour term of about **+0.027 in G
and −0.026 in R mag/mag**, equal and opposite. That is a model failure, not a measurement.

Two things were ruled out before changing the model:

- **The edge width.** The 10–90% width of 15 nm was assumed from the S50, and a shallower
  edge extends the filter's own red tail. Scanned from 5 to 110 nm, the G/R gap stays
  between +24 and +27 nm and widens slightly. Not it.
- **The pipeline.** Two unrelated sensors, one Sony and one OmniVision, producing
  near-identical systematics is not what two independent digitisation errors look like.
  See the control below — it exonerates the pipeline, which is what makes the gap real.

## The fix: fit the edge where the data is strong

All three planes see the IR cut, but only the **red** plane's *peak* sits under it — 70 to
85% QE right where the filter turns over. Its crossing is set by the edge and almost
nothing else. Green and blue are in a deep trough there (14–33% and 2–5%), which is the
thinnest and steepest part of a printed curve and the easiest few percent to digitise
wrong.

So the edge is measured from R alone, and the residual is attributed to the weak shoulders
of G and B as a multiplicative scale, ramped in past each dye's own peak so it never
touches the band's core. This is the same correction the S50 already needed for blue
(×1.15), applied to green as well.

| | IR-cut red edge | G shoulder | B shoulder |
|---|---|---|---|
| S50 *(control, not shipped)* | 692.8 nm | ×1.016 | ×1.204 |
| **S30 Pro** (IMX585) | **678.4 nm** | **×1.723** | **×2.452** |
| **S50 Pro** (OS08B10) | **684.9 nm** | **×1.644** | **×2.043** |

Blue edge 390 nm, 10–90% width 15 nm, peak 0.97, as for the S50. After fitting, every
residual colour term on every plane of every unit is zero to ±0.001 mag/mag.

```{figure} figures/unit-passbands.png
:alt: Total-system passbands for the S50, S30 Pro and S50 Pro, with each fitted IR-cut edge marked.

The three models' total-system passbands, with each unit's fitted IR-cut red edge. The
two newer models cut about 10 nm bluer than the S50 and have noticeably more green and
blue response under the red plane.
```

## The control

Everything above rests on the pipeline being trustworthy, so it was run on the S50 first —
the same code, the same offline catalogue, the same synthesised colour, on the 15-minute
local stacks the published S50 fit used. The published values came from a different script
in a different repository, against a TAP catalogue and real Gaia colours.

| | this pipeline | published |
|---|---|---|
| G crossing | 693.3 nm | 692 nm |
| R crossing | 692.8 nm | 693 nm |
| adopted edge | 693.1 nm | 693.0 nm |
| B shoulder | ×1.204 | ×1.15 (a three-stack scan gave 1.18) |
| **G shoulder** | **×1.016** | *(no such parameter)* |

The last row is the important one. The added green parameter comes out at **1.016 on the
S50** — inert exactly where the old model already worked, and departing from 1 only on the
sensors where it was broken. That is the signature of a parameter the data wanted, rather
than one bolted on to absorb a discrepancy.

## What it buys

Calibrator scatter per plane, against the nearest classical JKC band and against each
unit's own fitted band. Same calibrators, same frames, 3σ-clipped, per-frame median zero
point.

| unit | plane | JKC | own band | |
|---|---|---|---|---|
| S50 *(control)* | B | 360 mmag | **51** | 7.1× |
| | G | 64 | 40 | 1.6× |
| | R | 94 | 44 | 2.1× |
| **S30 Pro** | B | 263 | **81** | 3.2× |
| | G | 58 | 56 | 1.04× |
| | R | 131 | 83 | 1.6× |
| **S50 Pro** | B | 317 | **89** | 3.6× |
| | G | 87 | 76 | 1.15× |
| | R | 115 | 82 | 1.4× |

The pattern matches the S50's: blue gains most (3–4×), red gains usefully, green barely —
green was already close to Johnson V, which is why it is the science band. The absolute
floors are higher than the S50's because these are shorter stacks (260–380 s of on-sky
time against 900 s), not because the passbands are worse.

## What to be sceptical of

**The shoulder corrections are large.** ×1.7 in green and ×2.0–2.5 in blue means the
digitised curves understate those channels' red response by about a factor of two. On the
S50 the same correction was ×1.15 and ×1.02. The fits are well-determined and the scatter
improves, but "the datasheet understates a thin printed line" does not stretch to a factor
of two.

A mechanism consistent with the data, and worth testing rather than believing: **optical
and electrical crosstalk between neighbouring pixels**. A long-wavelength photon is
absorbed deep in the silicon and the carrier can be collected by a neighbouring pixel, so
the green and blue channels acquire apparent red response that no dye transmission curve
shows. It grows with wavelength and with smaller pixels, and both new units are 4K sensors
where the S50 is 1080p. That would raise G and B in exactly the 600–800 nm window the fit
is pushing up, and would leave R — where the signal is the dye's own — alone.

Until that is measured, treat the shoulder scales as *effective* parameters that make the
photometry work, not as a claim about the physical filter. They are fitted on one field at
one airmass range, and the two S30 Pro units and three S50 Pro units were not fitted
separately, so a unit-to-unit spread is not yet known.

**Nothing here is in the catalogue.** `gaiadb` carries `r/g/b_s50_mag` only. Adding
`*_s30pro` and `*_s50pro` columns means another schema version and another pass over the
114 GB of bulk XP spectra; `tools/build_gaia_s50.py` is the tool that would do it.
