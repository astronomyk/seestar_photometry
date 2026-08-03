# Light curves

## The two-table data model

`pipeline.build_measurements` writes two tables, and the split matters:

**`stars.ecsv`** — one row per catalogue source measured: position, magnitudes, colour,
variability flag, separation from the target, `is_target`. Small. This is the table you
*query* to choose a comparison ensemble.

**`measurements.ecsv`** — the long/tidy table, one row per `(source_id, frame, band)`:
forced-aperture flux, instrumental and calibrated magnitude, mid-exposure MJD/BJD, airmass,
SEP flag, `on_chip`, `max_pix_value`, and the frame's zero point. Large, and never queried
directly for selection — only reduced.

`frames.ecsv` joins to `measurements` on `frame` when you need per-frame quality.

## Forced photometry, not detect-then-match

Every catalogue position is measured in every frame, whether or not a source was detected
there. Consequences:

- a star keeps a row in **every** frame it lands on, so series never go ragged;
- the comparison ensemble is **identical** frame to frame.

With detect-then-match, a faint comparison drops out exactly on the frames where conditions
were poor — precisely when you need it — so the ensemble composition changes underneath the
light curve, and the change correlates with the conditions you are trying to divide out.

A consequence to expect rather than fix: **non-positive fluxes are normal**. A faint or
absent source measured in a fixed aperture on a background-subtracted image scatters either
side of zero. `photometry.instrumental_mag` turns those into `nan` rather than raising.

## Timing

Getting mid-exposure right matters as much as the photometry for anything with structure on
tens of minutes — a stack is a 10–15 minute integration, so a mid-point error of half that
smears and shifts a transit.

`lightcurves.frame_times` prefers, in order:

1. **The true exposure span**, `OB-START`/`OB-END` (CrowdSky). Exact, including inter-sub
   overhead. Recorded as `time_source="span"`.
2. **`DATE-OBS + total_exptime/2`** (native Seestar, which carries no end time). This ignores
   overhead between subs — unknowable from the header — so it runs slightly **early**: for a
   390 s on-sky stack spanning ~650 s of wall clock, by ~2 minutes.

`bjd_tdb` is the barycentric Julian date at mid-exposure. The correction uses the *target*
position: over a sub-degree field the per-star light-travel-time difference is < 1 s, so one
coordinate serves the whole frame. It comes back `nan` if the header has no site.

IERS auto-download is disabled and a stale bundled table is extrapolated without complaint —
far more accuracy than any of this needs.

## Choosing a comparison ensemble

`lightcurves.select_comparisons` filters on what actually matters:

- **proximity** — shared atmosphere, and on an Alt-Az mount a shared amount of field
  rotation, so the position-dependent aperture loss largely cancels;
- **similar brightness** — similar noise, away from both saturation and the detection floor;
- **similar colour** — minimises differential-extinction and colour-term residuals;

and drops catalogue-flagged variables.

A larger ensemble is **not** automatically better: adding faint or distant stars adds noise
and systematics faster than it averages them down. Somewhere around 10–30 well-matched
comparisons is usually the sweet spot.

Prefer setting `Target.source_id` explicitly over letting `build_stars` flag "the source
nearest the pointing". That heuristic is right for a field centred on an isolated target and
wrong in a crowded one.

## The ensemble differential

Each frame gets an **ensemble zero point**: the mean over its valid comparisons of
`(catalogue mag − instrumental mag)`. The target's calibrated magnitude is then
`m_target + ensemble ZP`. This removes the per-frame transparency (common to target and
comparisons, so it cancels) *and* ties the target to the comparisons' catalogue scale in one
step.

**Referencing each comparison to its own catalogue magnitude is the load-bearing detail.**
It is not a mean flux ratio. Because of it:

- a comparison dropping in or out between frames does **not** shift the zero point;
- the scatter of the per-comparison zero points is the genuine measurement error, not the
  intrinsic brightness spread of the ensemble.

`tests/test_lightcurves.py::test_comparison_dropout_does_not_shift_the_zeropoint` asserts
exactly this: a comparison vanishing halfway through must not step the curve. With a
mean-flux-ratio ensemble it would, and the light curve would show a discontinuity precisely
at the dropout.

Only valid measurements are used: on chip, finite magnitude, SEP `flag ≤ max_flag`. Frames
with fewer than `min_comp` (default 3) valid comparisons are dropped — below that the
ensemble is not averaging anything and its error is unmeasurable.

## Check the comparisons before believing the target

`lightcurves.comparison_curves` builds a differential curve for **each comparison, treated
as the target** and measured against all the others. This is the single most informative
check available, and it is worth reading *before* the target's curve.

A genuinely constant star gives a flat curve at the noise floor. A variable, a blend, a
saturated star or one drifting off the chip stands out immediately — and left in, it inflates
the target's error bars or, worse, imprints its own variability on the target. The
`lc_comparison_grid` diagnostic sorts panels worst-first and marks anything above 50 mmag
with a reserved status colour *and* the word CHECK.

`lc_noise_floor` plots achieved scatter against comparison brightness: bright stars sit on a
systematic floor, faint ones climb along photon noise. Where the floor sits **is** the
precision of the dataset.

## The MW Cam result

The c17 Alt-Az dataset recovers MW Cam's δ-Scuti pulsation at **P = 0.1294 d**
(Lomb-Scargle, FAP ~3e-11) with a scatter floor of **~23 mmag**.

Two findings about that floor:

- An aperture sweep put its minimum at **0.95 enclosed flux**, not 0.90 — a larger aperture
  damps the position-dependent aperture loss from **Alt-Az field rotation**. Hence the
  separate `enclosed_lightcurve` knob.
- Comparing the same target in EQ mode against Alt-Az showed a **~0.4 mag depth penalty**
  from field-rotation PSF broadening in Alt-Az. Mount mode is worth recording (`eqmode`) and
  worth checking before comparing datasets.

## Period analysis

`lightcurves.periodogram` (Lomb-Scargle) and `phase_fold`. Watch for aliases: a
nightly-cadence dataset has strong power at 1 day and its harmonics, so a period close to an
integer fraction of a day needs the phase fold before it is believed. The `lc_phase_fold`
panel draws two cycles with a binned mean so the shape reads continuously across the wrap.

## Targets on extended emission

For a supernova in a galaxy, or a star on a nebula, `contamination` estimates the host flux
inside the aperture **empirically**: sample the same aperture at the target's galactocentric
radius over clean azimuths and take the median.

Empirical rather than model-based because real hosts have structure at a given radius —
NGC 3310's circumnuclear starburst ring sits almost exactly at the radius of SN 2026sqf, and
a smooth Sérsic model under-predicts it substantially.

Two assumptions, both reported back so you can tell when they fail: the host is roughly
axisymmetric at that radius, and enough azimuths are clean of stars. Check `n_azimuth` —
below `min_azimuths` the correction degrades to zero rather than to a wild value. Note the
star mask deliberately **keeps** the galaxy nucleus, since masking it would remove the
emission being measured.
