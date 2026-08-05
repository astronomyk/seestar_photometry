# Calibrating a frame: the zero point

Turning aperture counts into magnitudes needs one number per frame, and it has to be
re-derived per frame because transparency changes. The package fits

$$V = m_\mathrm{inst} + \mathrm{ZP} + k\,(B-R - \mathrm{colour_0})$$

against Gaia DR3 **synthetic** Johnson-Kron-Cousins V.

## Do it

```python
from seestar_photometry import astrometry, calibration, examples, photometry

frame = examples.stack()
ext = photometry.extract_sources(frame)
ext.match_gaia(examples.gaia(), wcs=examples.wcs())

cal = calibration.fit_zeropoint(ext.sources, band="G")
print(f"ZP   = {cal.zeropoint:.3f} +- {cal.zeropoint_err:.3f}")
print(f"k    = {cal.colour_term:+.3f} at colour0 = {cal.colour0:.3f} ({cal.colour_label})")
print(f"rms  = {cal.rms:.3f} mag over {cal.n_stars} stars")
```

```
ZP   = 23.644 +- 0.003
k    = +0.036 at colour0 = 1.137 (B-R (JKC))
rms  = 0.016 mag over 59 stars
```

## Why each piece is the way it is

**Gaia synthetic V, not APASS.** It is homogeneous across the whole sky, and it carries
colour information, which is what makes a per-frame colour term fittable at all. Sources
without synthetic photometry come through masked, so a non-masked `v_jkc_mag` always means a
usable calibrator.

**A colour term, because green is not V.** The green Bayer channel through the Seestar's
IRCUT filter is a good but inexact match to Johnson V. The residual mismatch is small
(k = +0.036 here) but it is real, and it is colour-dependent, so ignoring it puts a
colour-correlated error into every magnitude.

**ZP quoted at `colour0`, the field's median colour.** This decorrelates the zero point from
the colour term. Without it the two trade off against each other and per-frame zero points
are not comparable — which defeats the purpose of measuring one per frame.

**A magnitude window, V ∈ [10, 14].** Above it stars saturate and the flux stops tracking V;
below it the noise floor fans the relation out. Either intruding into the fit biases it.

**Sigma clipping, and variables excluded.** Catalogue-flagged variables are dropped, and
3σ outliers are clipped over two refits.

## Reading the result

| Field | Here | What it tells you |
|---|---|---|
| `rms` | 0.016 | The master quality number. < 0.06 is photometric-grade; 0.1–0.3 marginal; > 0.5 junk. |
| `n_stars` | 59 | Calibrators surviving every cut. If this collapses, distrust the frame whatever its `rms` — with few stars the clipping has nothing to work with. |
| `chi2_red` | 214 | Scatter against the *photon-noise* prediction. **Not** expected to be ~1: Seestar frames sit on a ~0.03 mag systematic floor. Use it comparatively. |
| `zeropoint_err` | 0.003 | Formal error on the mean. Optimistic — `rms` is the honest per-star scatter. |

## Applying it

```python
mag = calibration.apply_calibration(m_inst, cal, colour=star_bp_rp)
```

With no colour for a source the colour term is dropped, which is exact at the field's median
colour and degrades linearly away from it — fine for a red-ish target, worth a thought for a
very blue one.

## Graceful degradation

The colour source cascades **B−R (JKC) → BP−RP → none**, and with no colour at all the fit
becomes a sigma-clipped zero-point mean with `colour_term = 0`. The reference magnitude
falls back from `v_jkc_mag` to a plain `v_mag`. That exists so the fit keeps working against
a different catalogue schema, not to squeeze out accuracy.

## Check it by looking

```python
from seestar_photometry import plots

plots.zeropoint_vs_colour(cal)          # the fit itself
plots.reference_vs_instrumental(cal)    # is the relation linear over the window?
plots.residual_vs(cal, against="v")     # any brightness-dependent bias?
plots.residual_map(cal, shape=frame.shape)   # spatial systematics
```

The residual map is the one people skip and shouldn't. Coherent structure across the frame —
a gradient, a corner, a ring — will **not** cancel differentially unless your comparison
stars happen to sit where your target does. It is what justifies a proximity cut when
choosing an ensemble.

Two stacks from the same telescope on the same night agree closely, as they must — the zero
point is a property of the unit and the sky, not of the integration:

```python
shallow = ...   # 760 s  -> ZP 23.644
deep    = ...   # 1460 s -> ZP 23.656
```

0.012 mag apart. A drift here would mean something has become exposure-dependent that
should not be.

More detail in [astrometry-and-gaia](../astrometry-and-gaia.md).
