# A differential light curve

The headline capability: ~20 mmag precision on a bright target, which is enough for δ-Scuti
pulsations, supernovae and exoplanet transits.

The example data has only three epochs, so what follows demonstrates the *mechanics*
faithfully but cannot demonstrate a real signal. Numbers from the full MW Cam dataset are
quoted where they matter.

## The two tables

```python
from seestar_photometry import examples, lightcurves, pipeline

proj = examples.project("work/")
pipeline.solve_all(proj)                       # no-op: solved sidecars come with the data
stars, meas = pipeline.build_measurements(proj)
print(len(stars), "sources,", len(meas), "measurements")
# 621 sources, 3804 measurements
```

`stars` is one row per source — position, magnitudes, colour, variability flag, separation
from the target. Small; this is the table you *query*.

`meas` is one row per `(source, frame, band)` — flux, instrumental and calibrated magnitude,
mid-exposure BJD, airmass, flags. Large; you reduce it, never query it for selection.

## Forced photometry, not detect-then-match

Every catalogue position is measured in every frame it lands on, whether or not anything was
detected there. That matters more than it sounds:

- a star keeps a row in **every** frame, so series never go ragged;
- the comparison ensemble is **identical** frame to frame.

With detect-then-match, a faint comparison drops out exactly on the frames where conditions
were poor — precisely when you need it — so the ensemble changes underneath the light curve,
correlated with the very conditions you are trying to divide out.

A consequence to expect rather than fix: **non-positive fluxes are normal**. A faint or absent
source in a fixed aperture on a background-subtracted image scatters either side of zero.
Those become `nan` magnitudes, not errors.

## Choosing comparisons

```python
comps = lightcurves.select_comparisons(
    stars,
    dmag=1.0,               # within +-1 mag of the target
    colour_tol=0.3,         # within +-0.3 in BP-RP
    max_sep_arcmin=15,      # close on the sky
)
```

The three cuts are not arbitrary:

**Proximity** — a shared atmosphere, and on an alt-az mount a shared amount of field
rotation, so the position-dependent aperture loss largely cancels.
**Similar brightness** — similar noise, and away from both saturation and the detection floor.
**Similar colour** — minimises differential-extinction and colour-term residuals.

Catalogue-flagged variables are dropped.

```{admonition} A bigger ensemble is not automatically better
:class: note
Adding faint or distant stars adds noise and systematics faster than it averages them down.
10–30 well-matched comparisons is usually the sweet spot. And a bright target is awkward: MW
Cam at V = 9.2 has almost nothing within ±1 mag, so a wider `mag_range` is unavoidable — at
which point the per-comparison check below stops being optional.
```

## The ensemble differential

```python
target_id = lightcurves.target_id_of(stars)
lc = lightcurves.differential_lightcurve(meas, target_id, comps, band="G")
print(lc["time", "mag", "mag_err", "n_comp"])
```

```
     time         mag    mag_err  n_comp
2461163.45220   9.2518   0.0053      62
2461163.48730   9.3189   0.0057      55
2461218.45136   9.3041   0.0091      57
```

Each frame gets an **ensemble zero point**: the mean over its valid comparisons of
`(catalogue mag − instrumental mag)`. The target's magnitude is then
`m_target + ensemble ZP`. This removes per-frame transparency and ties the target to the
comparisons' catalogue scale in one step.

**Each comparison is referenced to its own catalogue magnitude** — not to an ensemble mean
flux. That is the load-bearing detail. Because of it:

- a comparison dropping in or out between frames does **not** step the zero point;
- the scatter of the per-comparison zero points is the genuine measurement error, rather
  than the intrinsic brightness spread of the ensemble.

Note `n_comp` varies (62, 55, 57) as sources fall in and out of each frame's footprint. With
a mean-flux-ratio ensemble that variation would imprint itself on the curve; here it does not.
The package tests this explicitly by dropping a comparison halfway through a synthetic series
and asserting the curve does not step.

## Check the comparisons before believing the target

```python
curves = lightcurves.comparison_curves(meas, comps)
for sid, curve in sorted(curves.items(), key=lambda kv: -kv[1].meta["scatter"])[:3]:
    print(sid, f"{curve.meta['scatter']*1000:.0f} mmag")
```

Each comparison is measured against all the *others*. A constant star gives a flat curve at
the noise floor; a variable, a blend, a saturated star or one drifting off the chip stands
out immediately — and left in, it inflates the target's error bars or imprints its own
variability on the target.

On the three example epochs the spread runs from 0 to 68 mmag, but with three points that is
mostly meaningless — a two-point-plus-one series can be exactly flat by accident. **This
diagnostic needs tens of epochs to be worth reading.**

```python
from seestar_photometry import report
report.lightcurve_report(lc, stars, meas, comps, "diagnostics/")
```

`lc_comparison_grid.png` is the one to open first: panels sorted worst-first, with anything
above 50 mmag flagged in a status colour and labelled CHECK.

## Timing

```python
print(meas["bjd_tdb"][0])
```

Mid-exposure barycentric Julian date, TDB. This matters as much as the photometry for
anything with structure on tens of minutes: a stack *is* a 10–15 minute integration, so a
mid-point error of half that smears and shifts a transit.

The package prefers the true exposure span (`OB-START`/`OB-END`, written by CrowdSky) and
falls back to `DATE-OBS + total_exptime/2` for native frames, which runs ~2 minutes early
because inter-sub overhead is unknowable from the header. `lc.meta`/`time_source` records
which path was taken.

## Periods

```python
pg = lightcurves.periodogram(lc, min_period=0.02, max_period=1.0)
phase, mag = lightcurves.phase_fold(lc, pg["best_period"])
```

```{admonition} Check aliases before believing a period
:class: warning
On the full MW Cam c17 set — 52 frames over 1.32 days — the top peak came out at 0.11457 d,
while the true period is 0.1294 d. The two are exactly **1.0003 cycles/day apart**: a 1-day
alias, unresolvable at that baseline (frequency resolution 0.758 c/d). The real period was
present at power 0.246 against the alias's 0.291.

Always look at the phase fold. A real period folds into a coherent, repeatable shape; an
alias folds into scatter.
```

For reference, the full c17 dataset recovers MW Cam's δ-Scuti pulsation at **P = 0.1294 d**
with a scatter floor of about **23 mmag** — and reaching that floor took a tuned comparison
ensemble, not the first one selected. An aperture sweep also put the floor's minimum at
**0.95** enclosed flux rather than 0.90, because a slightly larger aperture damps the
position-dependent aperture loss caused by alt-az field rotation. That is why
`Project.enclosed_lightcurve` defaults to 0.95 while frame characterisation uses 0.90.

## Targets on extended emission

For a supernova in a galaxy, host light inside the aperture must be subtracted:

```python
from seestar_photometry import contamination

mask = contamination.star_mask(data_sub, rms, nucleus_xy)
host = contamination.galaxy_contamination(data_sub, nucleus_xy, sn_xy, aperture, mask=mask)
flux_corrected = flux - host["adu"]
```

The estimate is **empirical** — the same aperture sampled at the target's galactocentric
radius over clean azimuths — because real hosts have structure at a given radius. NGC 3310's
circumnuclear starburst ring sits almost exactly at the radius of SN 2026sqf, and a smooth
Sérsic model under-predicts it substantially. Check `host["n_azimuth"]`: below the minimum the
correction degrades to zero rather than to a wild value.

Full detail in [light-curves](../light-curves.md).
