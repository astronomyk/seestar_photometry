# seestar-photometry

Time-domain photometry from ZWO Seestar smart-telescope frames, calibrated onto Gaia DR3
synthetic Johnson V.

A Seestar is a $500 alt-azimuth telescope with a 50 mm aperture and a Bayer sensor. It is
not built for photometry. Used carefully it nonetheless delivers **~20 mmag differential
precision** on a bright target and reaches **V ≈ 17.5 in 15 minutes** — enough for
δ-Scuti pulsations, supernova light curves and exoplanet transits.

This package is the careful part. It grew out of two working analyses — a δ-Scuti light
curve of MW Cam and a Type II supernova in NGC 3310 — and encodes what those needed:

- **All three FITS layouts** the ecosystem produces: raw Bayer subs, on-board stacks, and
  CrowdSky multi-extension frames. They load to one array shape, and the two conflicting
  header dialects are resolved for you.
- **A per-frame, per-band aperture**, because the Seestar's PSF is chromatic and a shared
  radius biases every colour.
- **A re-solved WCS**, because the on-board one is off by ~1 arcminute and useless for
  photometry.
- **Ensemble differential photometry** that stays honest when comparison stars come and go.
- **Diagnostic figures for every stage**, so you can see that it worked rather than hope.

```{admonition} Every example here is executable
:class: tip
Every code block runs on real Seestar frames — cutouts, solved WCS sidecars and a trimmed
Gaia table — so no plate solver and no catalogue query is needed. `examples.download()`
fetches them once (~18 MB) and caches them; the wheel itself stays under 100 kB.
```

## Thirty seconds

```python
from seestar_photometry import examples, photometry, calibration

frame = examples.stack()                          # a real 760 s Seestar stack
ext = photometry.extract_sources(frame)           # detect + measure, per band
ext.match_gaia(examples.gaia(), wcs=examples.wcs())
cal = calibration.fit_zeropoint(ext.sources, band="G")

print(f"{len(ext.band('G'))} stars, ZP = {cal.zeropoint:.3f} ± {cal.rms:.3f} mag")
# 362 stars, ZP = 23.644 ± 0.016 mag
```

## Where to go

```{toctree}
:maxdepth: 1
:caption: Getting started

install
quickstart
```

```{toctree}
:maxdepth: 1
:caption: Use cases

usecases/stacking
usecases/zeropoint
usecases/saturation
usecases/depth
usecases/lightcurve
usecases/batch
```

```{toctree}
:maxdepth: 1
:caption: Reference

api
example-data
```

```{toctree}
:maxdepth: 1
:caption: Why it works this way

photometry-design
data-format
astrometry-and-gaia
frame-table
light-curves
diagnostics
architecture
migration-from-mwcam
```

The last section is the set of decision records: *why* each numerical choice is what it is,
with the measurements behind it. Read `photometry-design` before changing an aperture.

## Citing and caveats

Two numbers in particular are easy to over-trust, and both are documented in
[frame-table](frame-table.md):

- `v_lim_5sigma` runs **~0.7 mag optimistic** in absolute terms, because the aperture-noise
  model assumes independent pixels and demosaicing correlates them. It is a good *relative*
  metric within one processing chain and should not be compared across chains.
- `chi2_red` sits at 100–200 even on excellent frames. That is a real ~0.03 mag systematic
  floor, not a broken fit.
