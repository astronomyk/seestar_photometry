# The Seestar S50 photometric system

The S50 is a colour camera behind a fixed UV/IR-cut filter. Its three Bayer planes are a
photometric system in their own right — a system nobody has calibrated, because nobody
sells standard stars in it. This page measures where that system sits against
Johnson-Kron-Cousins and SDSS, gives the conversions between them, asks whether those
conversions need more than a linear colour term, and shows what calibrating against the
S50's own bands does to the zero point of a real stack.

Everything here comes from **one source and one integrator**. Every magnitude — S50, JKC
and SDSS alike — is an AB integral of the same Gaia DR3 XP spectrum through a different
passband:

    m = -2.5 log10[ INT f_lam S(lam) lam dlam  /  INT f_AB(lam) S(lam) lam dlam ]

so no cross-survey zero point, aperture correction or standardisation difference can leak
into a comparison. The sample is **235,740 stars** from 24 bulk XP files spread evenly
through the HEALPix ordering, which spreads them over Galactic latitude and therefore over
reddening. Rebuild it with `tools/build_s50_system_figures.py`.

:::{admonition} How far to trust the synthesis
:class: note

The same machinery, pointed at Johnson V, reproduces Gaia's own standardised
`v_jkc_mag` with a scatter of **8.2 mmag** over 231,924 stars, after removing a constant
+0.0364 mag offset that is just the AB-to-Vega difference for V. That 8 mmag is the floor
on everything below: it is the disagreement between this integration and an independent,
carefully standardised one, and no transformation here is better than it.
:::

## Where the three bands sit

```{figure} figures/passbands.png
:alt: S50 R, G and B passbands drawn against Johnson-Kron-Cousins BVRI and SDSS gri.

The S50's total-system response (filled, coloured) against the classical and Sloan
systems (grey outlines). The S50 curves are the sensor's Bayer dyes times the modelled
UV/IR-cut; the red edge at 693 nm is fitted from our own frames, not taken from a
datasheet.
```

| band | λ_eff [nm] | FWHM [nm] | | band | λ_eff [nm] | FWHM [nm] |
|---|---|---|---|---|---|---|
| **B_S50** | 506.0 | 114 | | B_JKC | 446.0 | 96 |
| **G_S50** | 564.4 | 124 | | V_JKC | 556.4 | 84 |
| **R_S50** | 624.6 | 114 | | R_JKC | 650.7 | 150 |
| | | | | I_JKC | 790.4 | 106 |
| | | | | g_SDSS | 475.1 | 124 |
| | | | | r_SDSS | 620.4 | 114 |
| | | | | i_SDSS | 751.9 | 122 |

λ_eff is photon-weighted — `INT S lam^2 dlam / INT S lam dlam` — so it is the wavelength
the band actually observes at rather than the midpoint of its edges.

Three things follow, and they drive everything else on this page.

**The three bands are bunched.** B_S50 to R_S50 spans 506–625 nm, 119 nm end to end. The
SDSS g-to-r baseline alone is 145 nm, and B-to-R in JKC is 205 nm. The S50 measures colour
over a shorter lever arm than either standard system, so its colours are intrinsically less
sensitive — and any conversion *out* of it is an extrapolation in wavelength.

**G_S50 is nearly V, and R_S50 is nearly Sloan r.** G_S50 sits 8 nm from Johnson V, and
R_S50 sits **4 nm** from Sloan r. That second one was a surprise: the natural reading of
"the red plane" is Cousins R, but R_S50 is 26 nm blue of R_JKC and essentially on top of
r_SDSS. The fits below confirm it — R_S50 → r is the tightest conversion on the page.

**B_S50 is not B, and nothing is I.** B_S50 is 60 nm red of Johnson B and 31 nm red of
Sloan g, because the sensor's blue dye has a broad red tail that the IR-cut does not
remove until 693 nm. And both I_JKC (790 nm) and i_SDSS (752 nm) lie entirely beyond that
cut: the S50 has no counterpart to them at all, which is visible in the figure as two grey
curves with nothing underneath.

## Converting to JKC and SDSS

The transformations are expressed in **B_S50 − R_S50**, the S50's own colour. That is the
whole point of a transformation: a relation keyed on BP−RP or B−V would require the
observer to already have the photometry they are trying to obtain. Written as

    target - S50_band = c0 + c1 x + c2 x^2 + c3 x^3,     x = B_S50 - R_S50

```{figure} figures/transformations.png
:alt: Six density panels of standard-minus-S50 magnitude against S50 colour, with linear and cubic fits.

Each conversion over 235,740 stars, with the linear (dashed) and cubic (solid) fits. The
quoted rms is the sigma-clipped residual scatter. Note the vertical scales: the whole
`r − R_S50` panel spans 0.2 mag, while `B − B_S50` spans 1.8.
```

Linear coefficients, valid over **0.0 ≤ x ≤ 1.65** (the 1st to 99.5th percentile of the
sample):

| conversion | c₀ | c₁ | rms |
|---|---|---|---|
| B_JKC − B_S50 | −0.1290 | +0.8089 | 48.4 mmag |
| V_JKC − G_S50 | −0.0446 | +0.0703 | 11.5 mmag |
| R_JKC − R_S50 | −0.0003 | −0.1482 | 8.5 mmag |
| g_SDSS − B_S50 | −0.1436 | +0.4143 | 44.2 mmag |
| r_SDSS − G_S50 | +0.0217 | −0.4502 | 7.6 mmag |
| r_SDSS − R_S50 | −0.0366 | +0.0230 | 9.0 mmag |

Cubic, over the same range:

| conversion | c₀ | c₁ | c₂ | c₃ | rms |
|---|---|---|---|---|---|
| B_JKC − B_S50 | +0.0017 | +0.5632 | −0.2171 | +0.3022 | 14.6 mmag |
| V_JKC − G_S50 | −0.0266 | +0.0716 | −0.1260 | +0.0973 | 5.0 mmag |
| R_JKC − R_S50 | +0.0098 | −0.2436 | +0.1931 | −0.1007 | 6.9 mmag |
| g_SDSS − B_S50 | −0.0367 | +0.2659 | −0.3284 | +0.3370 | 8.8 mmag |
| r_SDSS − G_S50 | +0.0225 | −0.4942 | +0.1138 | −0.0641 | 4.8 mmag |
| r_SDSS − R_S50 | −0.0135 | −0.0299 | −0.0101 | +0.0346 | 4.0 mmag |

Both tables, plus the quadratic fits, are in
[`s50-transformations.ecsv`](s50-transformations.ecsv) next to this page.

All magnitudes are **AB**. Landing on the Vega-based JKC scale that
`gaiadb`'s `v_jkc_mag` uses means one further constant per band; for V that constant is
the +0.0364 mag measured above.

:::{warning}
Do not use these outside 0 ≤ B_S50 − R_S50 ≤ 1.65. That range is the 1st to 99.5th
percentile of the sample, measured rather than chosen: above 1.65 there are almost no
stars to constrain a cubic, and a cubic with nothing holding its red end down does not
flatten out, it turns over. Fitting the same data to x = 2.2 — only 0.5% more stars —
already produced corrections a magnitude beyond anything the data supported. If your star
is redder than 1.65, the honest answer is that this page does not cover it.
:::

## Does it need more than a linear colour term?

```{figure} figures/polynomial-order.png
:alt: Residual rms against polynomial degree for each conversion, in two panels for JKC and SDSS.

Residual scatter against polynomial degree. The 10 mmag line is roughly where the
synthesis floor and typical per-frame photometric noise sit, so gains below it are
academic for a single frame.
```

**Yes for the blue plane, no for the red one.** The answer splits cleanly:

- **B_S50 → B_JKC or g_SDSS: a linear term is not enough.** Going linear → cubic takes
  B_JKC from 48.4 to 14.6 mmag and g_SDSS from 44.2 to 8.8 mmag, factors of 3.3 and 5.0.
  The relation is visibly curved in the figure above, and it curves for a physical reason:
  B_S50's red tail means the band's effective wavelength itself slides redward as a star
  gets redder, so the offset to a narrow blue band cannot be linear in colour.
- **R_S50 → r_SDSS: linear is already close, quadratic finishes it.** 9.0 → 4.2 mmag from
  degree 1 to 2, and degree 3 buys 0.2 mmag more. The bands nearly coincide, so there is
  little for the colour term to do.
- **R_S50 → R_JKC and G_S50 → V_JKC: modest gains, and the quadratic is not where they
  come from.** R_JKC goes 8.5 → 8.2 → 6.9 mmag; V_JKC goes 11.5 → 6.6 → 5.0. Note V's
  jump is at degree 2 while R's is at degree 3.

A practical reading: **use the cubic for anything involving B_S50, and the quadratic
elsewhere.** The extra terms cost nothing to evaluate, and the only reason to prefer the
linear form is if you are fitting the transformation yourself from a handful of stars.

### Where a polynomial stops helping

Look again at the `R − R_S50` panel in the transformation figure. Above x ≈ 1.05 the
relation **splits into two branches**, and the lower one falls away by 0.4 mag. The same
split appears in `V − G_S50`, and much more violently in the bands beyond the IR cut.

That is not fit error. It is M dwarfs separating from giants: below about 4000 K, TiO and
VO bands carve up the 600–800 nm region, and how much they bite depends on surface gravity
as well as temperature. A single colour cannot tell those two populations apart, so *no*
polynomial in B_S50 − R_S50 can follow both branches — the residual floor there is set by
the parameterisation, not by the degree.

The clearest evidence is the two bands with no S50 counterpart at all. Forcing a fit from
R_S50 across the IR cut gives:

| conversion | linear | quadratic | cubic |
|---|---|---|---|
| i_SDSS − R_S50 | 51.8 mmag | 42.8 mmag | 38.3 mmag |
| I_JKC − R_S50 | 64.6 mmag | 52.5 mmag | 47.0 mmag |

Three degrees of freedom buy a 27% improvement on a scatter that is an order of magnitude
worse than the in-range conversions. **The S50 cannot be transformed into I or i.** If you
need those bands, you need a different filter, not a better polynomial.

## What this does to a zero point

The reason any of this matters is calibration. The standard path fits

    V = m_inst + ZP + k (colour)

against a reference catalogue, and `k` exists to absorb the mismatch between the
instrument's band and the reference's. If the reference *is* the instrument's band, there
is nothing for `k` to absorb.

Here is one 15-minute 20 Cep stack, each Bayer plane calibrated twice against the same
calibrator stars — once against the closest classical band, once against the synthetic S50
band from `gaiadb`. The zero point is a median in both cases, so the residual slope is
what the reference leaves behind rather than something the fit has already removed.

```{figure} figures/frame-zeropoint.png
:alt: Zero-point residual against BP-RP for three Bayer planes, comparing classical and S50 references.

`n2_c17_dither` 15-minute stack, 238 calibrators, V ∈ [9.5, 13.0), SNR > 30, unblended,
unsaturated, non-variable. The colour axis is BP−RP — deliberately *not* an S50 colour,
which would be testing the S50 reference against itself. Scatter is 3σ-clipped.
```

| plane | classical reference | S50 reference | |
|---|---|---|---|
| **B** | 348 mmag, slope **+0.577** | **43 mmag**, slope +0.028 | 8.0× tighter |
| **G** | 57 mmag, slope +0.057 | 45 mmag, slope +0.001 | 1.3× tighter |
| **R** | 86 mmag, slope −0.126 | 42 mmag, slope +0.003 | 2.0× tighter |

The B plane is the headline: against Johnson B its residuals run up a **+0.58 mag/mag**
ramp across the field's colour range, and simply changing the reference collapses that to
+0.028 and cuts the scatter by a factor of eight. But the more telling number is the
right-hand column — **all three planes land on 42–45 mmag**, within 3 mmag of each other,
once each is calibrated against its own band. That common floor is the frame's own
photometric noise. Against the classical references the three planes are spread over a
factor of six, and the spread is entirely the mismatch between each Bayer plane and the
classical band it was being forced onto.

:::{admonition} The blend cut is doing a lot of this work
:class: caution

"Unblended" in the caption means no catalogue neighbour within 1.5 aperture radii that is
less than 2.5 mag fainter. Drop it and this frame's G-plane scatter goes from 45 to 100
mmag, and the deep stack below from 87 to 164. At 2.4 arcsec/px with a ~10 arcsec FWHM a
blended aperture is measuring two stars, and in a Galactic-plane field like 20 Cep it is
29% of the sources on the chip. It is the single easiest cut to leave out and the most
expensive.
:::

### Seen against magnitude instead of colour

The residual-against-colour view above answers "is the reference on the right system". The
complementary view — per-star zero point against reference magnitude — answers "does the
zero point mean the same thing at every brightness", and shows every way a calibration can
fail in one frame.

```{figure} figures/zp-linearity.png
:alt: Per-star zero point against reference magnitude for Gaia synthetic V and G_S50, with a residual panel.

A 60-minute 20 Cep stack (141 subs, 2820 s), green plane, calibrated against Gaia
synthetic V and against `g_s50_mag`. Shaded band is the calibrator window; the lower
panel is the calibrators about each reference's own zero point. A deeper stack than the
one above, on purpose: it reaches ~2.5 mag fainter, which is where the two references
part company.
```

Three regimes, left to right:

- **Brighter than V ≈ 9.5 the zero point falls away** by up to half a magnitude. That is
  flux clipped in the individual subs. Note that no ADU test catches it: this stack's
  brightest pixel is 53,140 ADU, comfortably under any sensible fraction of the 65,535
  clip level, because averaging 141 subs pulls the clipped peak down while the lost flux
  stays lost. **On a stack, saturation has to be excluded by magnitude, not by pixel
  value** — which is what the calibrator window's bright edge is for.
- **Across the calibrator window the cloud is flat**, and this is where the two references
  separate: **52 mmag against synthetic V, 35 mmag against G_S50**, over the same 210
  stars. A third of the scatter that looked like photometric noise was the reference.
- **Fainter than V ≈ 13 the cloud flares upward and only upward.** The asymmetry is the
  giveaway that it is not symmetric measurement noise: a faint star's aperture can lose
  flux to the background subtraction and to a neighbour, but it cannot gain any, so the
  measured magnitude scatters faint and the per-star zero point scatters high.

One frame is an anecdote, so the same comparison over all 57 stacks:

```{figure} figures/survey-scatter.png
:alt: Strip plot of calibrator rms for 57 stacks, three planes, two references, with medians.

Every 15-minute stack of the three 20 Cep night-2 datasets, both references. Vertical
ticks are medians. The clumping within each group is the three datasets, which were taken
under different conditions.
```

| plane | classical (median) | S50 (median) | classical slope | S50 slope |
|---|---|---|---|---|
| B | 343 mmag | **84 mmag** | +0.541 | −0.010 |
| G | 91 mmag | 85 mmag | +0.054 | −0.002 |
| R | 126 mmag | 87 mmag | −0.131 | −0.003 |

All three planes land on the same **84–87 mmag** floor once they are calibrated against
their own band, and all three colour terms collapse to within 0.01 mag/mag of zero — the
same pattern as the single frame, over 57 of them.

These rms values are *unclipped* standard deviations, which is why they sit above the
3σ-clipped 42–45 mmag quoted for the single frame; the two columns are comparable to each
other but not to that table. What is comparable either way is the collapse of the spread
between planes, and the colour terms.

## Using this

`gaiadb` carries `r_s50_mag`, `g_s50_mag` and `b_s50_mag` for every catalogue source with
a Gaia sampled XP spectrum, so the S50-reference calibration above needs no transformation
at all — see [astrometry-and-gaia](../astrometry-and-gaia.md). The transformations on this
page are for the other direction: **reporting** an S50 measurement to someone who works in
V or r.

Three caveats worth carrying:

- The S50 columns exist only below **G = 15.0**, which is 23% of the catalogue. Fainter
  calibrators still have to go through `v_jkc_mag` and a fitted colour term.
- The passband model behind all of this has one fitted parameter that is not tightly
  pinned: the blue dye's red-tail scale, adopted at ×1.15 between a three-stack estimate of
  1.18 and a 57-stack estimate of 1.12. The B-plane transformations are the ones that
  depend on it.
- Nothing in the calibration path reads the S50 columns yet. `calibration.fit_zeropoint`
  still fits V with a colour term, which is why the numbers above were produced by
  `tools/build_s50_system_figures.py` rather than by the pipeline.
