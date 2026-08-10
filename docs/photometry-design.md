# Photometry design

How `photometry` measures flux, and why each choice was made. The numbers below were
measured on the June–July 2026 MW Cam datasets (S50 + S30pro).

## 1. SEP fixed circular aperture is the standard; Kron is a cross-check only

`extract_sources` detects with SEP (`sep.extract`, low threshold) and measures flux in a
**fixed circular aperture**, keeping only sources with **SNR > 5**. Kron/AUTO photometry
exists but lives in `kron`, and is not the science path.

Detection uses a deliberately **low threshold (2σ)** so the SNR > 5 cut on the *aperture*
flux — not the per-pixel detection step — defines the final sample.

**Why circle over Kron.** For differential/ensemble photometry of point sources, a
consistent aperture makes the aperture correction common to all stars and frames, so it
cancels in Δm and is absorbed into the per-frame zero point. Kron's per-source adaptive
aperture introduces source-dependent (and frame-dependent) systematics that don't cancel.

| Goal | Fixed circle | Kron |
|------|:---:|:---:|
| Differential light curves | ✅ apcorr cancels in Δm | ❌ per-source apcorr doesn't cancel |
| Detection limits / mag–SNR | ✅ single clean −0.4 locus | ❌ double-stream artefact |
| Gaia zero point | ✅ apcorr absorbed into ZP | ⚠️ fallback flux isn't truly total |
| Star-flat / ubercal | ✅ fixed aperture is the substrate | ❌ adaptive aperture hides spatial PSF |

**The Kron "secondary stream".** In a mag–SNR diagram Kron produces a tight secondary
streak offset from the main locus. Cause: sources whose Kron radius is too small to
integrate fall back to a **fixed minimum circular aperture** (1.75 px) — that
identical-aperture sub-population lines up on its own constant-area −0.4 line. `kron`
exposes a `kron_fallback` column; drop those rows before treating Kron flux as a total
magnitude. A fixed aperture collapses the whole diagram to one clean locus.

*Implementation note:* SEP's `sum_ellipse`/`kron_radius` require `|theta| < π/2` strictly,
but `extract` occasionally returns `theta` a float-epsilon outside that range.
`kron._kron_flux` clips theta just inside, otherwise it raises `invalid aperture
parameters` on a subset of frames.

## 2. Aperture sized per frame from the PSF, by enclosed flux fraction

The aperture radius is **not** hardcoded. It is the radius enclosing a given fraction of
the flux, read off a **per-frame** curve of growth (`curve_of_growth`, median over
bright/round/isolated stars). At 90% this lands at **~1.1–1.3 × FWHM**. Overrides:
`n_fwhm=` (a multiple of the measured FWHM) or an explicit `aperture=`.

**Why COG-driven, not a fixed multiple like 3.5 × FWHM.** Aperture sizing matters far
more than "it just shifts the zero point" suggests — that intuition holds for *bright*
stars (source-noise dominated) and for the per-frame ZP mean, but **not** for the
faint-end deliverables. Measured trade-off (S50, FWHM ≈ 4 px):

| Aperture choice | radius | in FWHM | faint-star SNR delivered |
|---|---|---|---|
| SNR-optimal | 2.5 px | 0.62 × FWHM | 1.00 |
| **90%-enclosed** | 4.7 px | 1.16 × FWHM | 0.80 |
| 3.5 × FWHM | 14 px | 3.5 × FWHM | 0.29 |

A 3.5 × FWHM aperture throws away ~70% of faint-star SNR, which would inflate the
reported **5σ detection limit by ~1.35 mag** — a headline number. So the detection-limit
work in particular must not use an oversized aperture. Enclosed-fraction sizing is the
compromise: near the SNR sweet spot while robust enough for differential work.

FWHM is derived as `2 × half-light radius` (`sep.flux_radius`, Gaussian assumption;
`GAUSS_FWHM = 2.0`).

### Two fractions, deliberately

`Project` carries **two** knobs, and they differ on purpose:

| Knob | Default | Used for | Why |
|---|---|---|---|
| `enclosed_characterise` | **0.90** | `frames.ecsv`, depth limits | Near the SNR sweet spot; the published depth numbers use it |
| `enclosed_lightcurve` | **0.95** | forced photometry for light curves | An aperture sweep found the differential scatter floor minimises here |

The light-curve value was raised from 0.90 after sweeping 0.85/0.90/0.95/0.99 against the
achieved scatter floor: a slightly larger aperture damps the *position-dependent* aperture
loss caused by Alt-Az field rotation, at negligible sky-noise cost for photon-dominated
targets. Do not "unify" these — they optimise different quantities, and
`tests/test_pipeline.py` asserts they stay distinct.

### A nebula will size your aperture if you let it

The COG sample is "bright, round, isolated". On a field containing extended emission those
three cuts do not select stars — they select *the nebula*. Measured on a real M27 frame
(600 s stack, 5469 green detections):

| | |
|---|---|
| bright + round (SNR > 50, b/a > 0.7) | 578 sources |
| median nearest-neighbour distance | **10.6 px**, against a 40 px isolation cut |
| sources clearing isolation | **1** — M27 itself, `a = 42 px`, `b/a = 0.82` |
| aperture from that one-object COG | **19.0 px** |
| aperture the stars actually wanted | **5.3 px** |

Two failure modes compound:

- **Isolation selects for extended objects.** In a rich field nothing *stellar* has 40 px of
  clear space, but a nebula does, precisely because it is large enough to have swallowed or
  outgrown its neighbours. Roundness does not save you: a planetary nebula is round.
- **A one-star COG returns a finite radius.** The documented "falls back to 1.2 × FWHM"
  only triggers on a *non-finite* result, so a garbage aperture measured from a single
  extended object was used silently, with `n_stars = 1` the only clue.

Both are now closed. `COG_MAX_SIZE_RATIO` (3.0) rejects candidates whose semi-major axis
exceeds 3× the median of the frame's bright round sources — a star sits at 1, M27 at 21 —
and `MIN_COG_STARS` (5) forces the FWHM fallback rather than trusting a COG built on
almost nothing. On top of that every entry point takes `mask=` (`True` = ignore, SEP's
convention) for emission you know about in advance; `photometry.sky_mask` builds one from
sky circles, and `Project.mask` is the per-frame callable hook the pipeline uses.

**The FWHM is not affected — only the aperture.** Masking M27 moved the measured green FWHM
by 0.01 px (4.41 → 4.40) while moving the 90% aperture from 19.03 px to 5.29 px. `measure_fwhm`
takes a median of per-source half-light radii, and a median of scalars shrugs off a minority
that a median of whole growth *curves* does not: one non-convergent curve depresses the
median curve at every radius, so the radius where it first reaches 0.90 slides outward.

Note what `mask=` deliberately does **not** do. It is not passed to `sep.extract`, because
that does not remove a bright extended source — SEP drops the masked pixels from the
footprint but the surrounding flux still clears the threshold, so one nebula fragments into
a ring of detections around the mask boundary, whose centroids all sit *outside* the mask.
Rejecting on the centroid of the un-masked segmentation removes the object in one piece.
And it never changes a reported flux: a star whose aperture overlaps the mask is still
summed over all its pixels, because silently altering a flux is worse than reporting a
contaminated one.

## 3. The aperture is sized per *band*, not pinned to green

Each of R, G, B is sized to **its own** enclosure radius. They are *not* all pinned to the
green radius.

**Why — the Seestar PSF is chromatic.** Measured FWHM ratios: `R/G ≈ 1.04–1.20`,
`B/G ≈ 1.06–1.30`. Both red **and blue** are broader than green. This is **chromatic
aberration**, not diffraction: the PSF (~9″) is ~10× the diffraction limit (~0.9″), so it
is seeing/optics/tracking-dominated, and the cheap optics focus best near green so both
sides defocus. (A diffraction-limited system would give a *narrower* blue PSF — the data
rules that out.)

If all bands shared the green radius each would enclose a different flux fraction → a
band-dependent aperture correction → a biased **B−R colour**, the quantity the Gaia
calibration depends on. Worse, the broadening tracks the auto-focus, so the bias **varied
frame to frame from +0.04 to −0.13 mag and changed sign** — it cannot be removed by a
single colour-term constant.

Sizing each band to the **same enclosed fraction** instead makes the aperture correction
identical across bands (−0.114 mag = `2.5·log₁₀(0.90)` by construction), so it cancels in
colours. This reversed an earlier, incorrect "a common radius keeps colours clean"
rationale: for a chromatic PSF it is a common *fraction*, not a common *radius*.

The `*_fwhm` panel and the `*_colour_fidelity` panel of the per-frame diagnostic set are
the two figures that check this is still true on a new dataset.

## 4. Aperture correction

`aperture_correction(frame, radius, cog=)` returns `2.5·log₁₀(frac(radius))` — the offset
from a fixed-aperture instrumental magnitude to a total magnitude (negative; total flux is
brighter). For relative work it is mostly a cross-check: it is absorbed into the
per-frame, per-band zero point. It is needed for absolute total magnitudes.

With per-band 90% sizing it is −0.114 mag in every band by construction. This absorption
is asserted directly in `tests/conftest.py::expected_zeropoint` — a fitted zero point comes
out exactly `2.5·log₁₀(enclosed)` below the true total-flux zero point.

## Open items and caveats

- **Per-frame apcorr genuinely varies** — from −0.10 to −0.69 mag at a *fixed* radius
  across frames of differing focus. The per-frame COG handles this; never reuse one
  frame's correction across a stack.
- **Intra-frame spatial PSF variation** is not addressed by aperture choice at all. That is
  the star-flat / ubercal problem. The `*_residual_map` diagnostic panel is what reveals
  it, and a proximity cut on the comparison ensemble is the current mitigation.
- The aperture is sized from **bright, round, isolated, point-like** stars. In sparse or
  badly focused frames the COG star count drops (often ~15–20) — still enough for a median,
  but worth watching via `cog.meta["n_stars"]`. Below `MIN_COG_STARS` (5) the code falls
  back to `1.2 × FWHM`.
- **In a crowded field the COG path may never engage at all.** The isolation requirement is
  `2 × COG_RADII.max()` = 40 px; on the M27 fields the median nearest-neighbour distance is
  ~10 px, so after the size cut correctly removes the nebula the sample can be empty and
  every band silently takes the `1.2 × FWHM` fallback — meaning `enclosed_lightcurve = 0.95`
  has no effect there. `curve_of_growth(..., isolation=)` and `Project.isolation` let you
  lower it; check `n_stars` to see whether you are getting a measured aperture or the
  fallback. Lower it knowingly: the outer radii are what normalise the curve.
