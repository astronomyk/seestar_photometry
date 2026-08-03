# Diagnostic figures

Every figure exists to answer a specific question. This is the catalogue: what each panel
shows, and what "wrong" looks like in it.

`plots` draws individual figures (`(data, ..., ax=None) -> ax`, no file I/O, composable into
a notebook). `report` composes them into named sets and saves PNGs to
`work_dir/diagnostics/`. Filenames are deterministic, so re-running overwrites rather than
accumulating and any figure can be asked for by name.

```python
pipeline.build_frame_table(proj, diagnostics=3)   # dataset panels + 3 per-frame sets
report.lightcurve_report(lc, stars, meas, comps, proj.diagnostics_dir)
```

`diagnostics=True` gives the dataset panels only; an `int` also saves that many per-frame
sets. Per-frame sets re-measure their frames in the calling process, so they work even when
every row was already cached — which is the usual case when you come back to ask why
something looks wrong.

## Per frame — `report.frame_report`

| File | Question | What wrong looks like |
|---|---|---|
| `*_cog` | Is the aperture sized sensibly? | A curve still climbing at the last radius — the "total" normalisation is contaminated by a neighbour, so the aperture was sized too large for the whole frame |
| `*_fwhm` | Is the PSF chromatic, as the per-band aperture assumes? | Ratios at 1.00 (per-band sizing buying nothing) or far above R/G ≈ 1.2, B/G ≈ 1.3 (bad focus, trailing) |
| `*_zp_relation` | Is V vs instrumental magnitude linear over the fit window? | Flattening at the bright end (saturation) or fanning at the faint end (noise floor) *inside* the shaded window — move the window |
| `*_zp_colour` | The colour-term fit itself | A slope that swings frame to frame — something varies that a constant colour term can't absorb, usually focus via the chromatic PSF |
| `*_residual_v` | Brightness-dependent bias? | Any slope: non-linearity, or saturation creeping into the window |
| `*_residual_snr` | Is the noise model right? | Structure means it isn't |
| `*_residual_radius` | Vignetting / flat-field residual? | A trend with radius — and on Alt-Az, field-rotation-dependent aperture loss |
| `*_residual_map` | Spatial systematics | Any coherent structure — a gradient, a corner, a ring. This will **not** cancel differentially unless comparisons sit where the target does; it is what justifies a proximity cut |
| `*_mag_snr` | Where is the SNR > 5 boundary, per band? | A second offset streak — two populations measured differently (the Kron fallback signature) |
| `*_background` | Did the polynomial absorb the sky gradient? | Large-scale structure left in the residual — the sky isn't quadratic, and the pedestal (hence sky brightness) is biased. Bright-star halos in the residual are fine, they belong to the sources |
| `*_detections` | Do astrometry, detection and cross-match agree? | Apertures off-centre from stars, or a systematic offset between apertures and catalogue rings — a bad WCS, far easier to see here than in any statistic |
| `*_match_sep` | How good is the solve? | A broad distribution filling the tolerance, or a peak pressed against it — the solve is wrong even though it "succeeded". Good is ~0.8″ median |
| `*_colour_fidelity` | Does the per-band aperture keep colours unbiased? | Slope away from 1, or curvature — the bands are *not* enclosing matched fractions, exactly the bias a shared radius introduces |

## Per dataset — `report.frames_report`

| File | Question |
|---|---|
| `frames_zp_vs_time` | Transparency trends and cloud, per unit. Units are *expected* to be offset from each other — read the shape within a series, not the spacing between them |
| `frames_rms_hist` | How much of the dataset is usable, against the 0.06 cut |
| `frames_chi2_hist` | Read the shape and pick out the far-right frames. Do **not** compare against 1 — see [frame-table.md](frame-table.md) |
| `frames_depth_vs_exptime` | Does depth follow √t? Flattening means a systematic floor took over |
| `frames_depth_vs_{sky_sb,fwhm_G,airmass}` | The condition drivers, with exposure divided out. Check the frames actually span a range of the driver before believing a fitted slope |
| `frames_{fwhm_G,sky_sb,airmass}_vs_time` | Conditions across the night — makes a feature in a light curve interpretable, or dismissable |
| `frames_coverage` | `n_cal` per frame. A collapse means a frame to distrust whatever its `rms` says |
| `frames_summary` | All of the above on one page |

## Per light curve — `report.lightcurve_report`

| File | Question |
|---|---|
| `lc_finder` | Is the ensemble what you meant to select? Are comparisons spread around the target rather than clustered on one side? |
| **`lc_comparison_grid`** | **Read this first.** Each comparison measured against the others. Every panel should be flat at the noise floor; a trend, a step or excess scatter is a comparison that doesn't belong. Sorted worst-first, with anything above 50 mmag marked in a status colour *and* labelled CHECK |
| `lc_noise_floor` | Achieved scatter vs comparison brightness. Where the floor sits **is** the dataset's precision. A bright star above it is a bad comparison; one climbing again at the very bright end is saturating |
| `lc_differential` | The result. Compare the quoted scatter against the median error bar: agreement means you're at the photon limit, scatter ≫ errors means a systematic remains |
| `lc_ensemble_zp` | The correction being applied, with the star-to-star spread. A frame where the spread balloons had a comparison misbehave, and the target's point inherits it |
| `lc_raw_flux` | The uncorrected target flux. Raw noise plus a clean differential curve is the expected outcome — and confirms the differential step is doing real work rather than flattening a genuine signal |
| `lc_periodogram` | The signal and its FAP. Check the peak against the sampling: nightly cadence puts power at 1 day and harmonics |
| `lc_phase_fold` | Confirmation. A real period folds into a coherent shape; an alias folds into scatter |
| `lc_host_cutout` | For a target on extended emission: the aperture and the azimuthal sampling ring. Few markers, or markers all on one side, means the contamination estimate rests on a poor sample |

## Styling, and why it looks like this

All in `_style`. Committed to a **light surface** — these end up in notebooks, papers and
PNGs, all of which are light. There is no dark variant to keep in sync.

**Band colours are domain-mandated.** R/G/B are drawn red/green/blue because here the colour
*is* the band name; relabelling them would mislead the reader. That is not free: red vs green
is the classic colour-vision-deficiency collision. The chosen steps (`#e34948` / `#008300` /
`#2a78d6`) measure a worst-pair CVD separation of **ΔE 7.2 (protan)** — inside the 6–8 floor
band, not clear of it. Alternatives were measured; every greener or lighter green traded
protan separation for a worse tritan one, so these are the best available under the
constraint.

The floor band is legal **only with secondary encoding**, so hue is never the sole carrier of
band identity: each band also has its own **marker shape** (`BAND_MARKER`) and **line style**
(`BAND_STYLE`), and multi-band panels are always legended. **Do not remove those when editing
a figure** — they are the accessibility mechanism, not decoration.
`tests/test_report.py::test_band_palette_uses_secondary_encoding` guards it.

Non-band series (units, datasets) use a validated categorical palette in fixed slot order —
assigned by position and never cycled, so a filter that drops a series doesn't repaint the
survivors — plus the same marker-shape channel, since a scatter with more than three series
cannot clear the all-pairs floors on hue alone.

Two further conventions worth knowing before editing:

- **Dense clouds are drawn small and translucent**, not as large opaque marks. A field is
  hundreds to thousands of stars — a density distribution, not a set of labelled series — so
  overplotting, not mark size, is the thing to solve. Legends bump `markerscale` to
  compensate.
- **Legends go outside, above the axes** on multi-series panels (`_style.legend(outside=True)`),
  which is the only reliably collision-free placement, and the panel title is re-padded to
  clear them. Call `legend()` *after* `set_title()`.
- Status colours (`good`/`warning`/`serious`/`critical`) are **reserved** — never reused as a
  series colour, and always paired with a text label so state is not carried by hue alone.

## Failure behaviour

A diagnostic set is most valuable exactly when something is wrong, which is also when an
individual panel is most likely to fail (an empty band, a missing column, a degenerate fit).
`report._panel` therefore catches a failing panel, draws the error text into that image, and
continues — one dead panel never costs you the other twenty. The tests check for these
placeholders by size, because "the files exist" would otherwise be a weak assertion.
