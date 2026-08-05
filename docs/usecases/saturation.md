# The bright limit: where saturation starts

Every frame has a magnitude above which flux stops tracking brightness, because the peak
pixel has hit the sensor ceiling. Measure it, or you will quietly put clipped stars into a
zero-point fit or use one as a comparison star.

## Do it

```python
from seestar_photometry import calibration, examples, photometry, quality

frame = examples.stack_saturated()      # a 280 s stack whose brightest star clips
ext = photometry.extract_sources(frame)
ext.match_gaia(examples.gaia(), wcs=examples.wcs("stack_saturated"))
cal = calibration.fit_zeropoint(ext.sources, band="G")

print(f"peak pixel  = {frame.g.max():.0f} ADU  (ceiling 65535)")
print(f"v_sat       = {quality.frame_quality(ext, cal)['v_sat']:.2f}")
```

```
peak pixel  = 64863 ADU  (ceiling 65535)
v_sat       = 8.09
```

Stars brighter than V ≈ 8.1 are unreliable in this frame.

## How it is measured, and why empirically

```python
v_sat = calibration.saturation_mag(ext.band("G"))
```

It is the **faintest** matched star whose peak pixel sits at or above the saturation level —
so everything brighter is suspect. Each source's peak is recorded during extraction as
`max_pix_value`, taken from the *raw* plane rather than the background-subtracted one, so it
can be compared against the ceiling directly.

The threshold is 95% of the 16-bit ceiling, not the ceiling itself:

```python
calibration.SAT_LEVEL     # 62258.25  ==  0.95 * 65535
```

A Seestar stack is a **scaled average** of its subs, so a pixel that hard-clipped in several
individual subs comes out below 65535 in the stack. Waiting for a true 65535 would miss most
real saturation; 0.95 catches the flat-topped stars in practice.

That is also why this is measured rather than derived. The clipping level in a stack depends
on how many subs contributed and how they were scaled, which is not something you can compute
from the header.

## `nan` is the right answer sometimes

```python
clean = examples.stack()          # nothing in this field clips
# ... extract, match, fit ...
quality.frame_quality(ext, cal)["v_sat"]
# nan
```

`nan` means no matched star in the field saturated — the bright limit is brighter than
anything present, so the frame cannot constrain it. Do **not** treat `nan` as "no saturation
problem" in general; it means "not measurable here". In the bundled saturated frame just
**1 of 266** green sources clips, which is all it takes to define the limit and all it takes
to poison a fit if you leave it in.

## Using it

The zero-point fit already protects itself: the V ∈ [10, 14] window sits well below any
plausible `v_sat`, so clipped stars are excluded by construction. Where you have to be
careful is choosing comparison stars for a light curve, and interpreting a bright target:

```python
frames_table = ...      # from pipeline.build_frame_table
usable = frames_table[frames_table["v_sat"] < target_mag]   # target fainter than the limit
```

```{admonition} A bright target is the awkward case
:class: warning
MW Cam is V = 9.2. On a long stack that is comfortably linear, but on a short one it can
approach `v_sat`. If your target is within ~1 mag of the frame's saturation limit, check
`max_pix_value` for the target itself rather than trusting the frame-level number, which is
derived from whichever star happened to be brightest.
```

## The other end

`v_lim_5sigma` is the faint limit, and `v_lim_100sigma` a high-SNR working limit. Together
with `v_sat` they bracket the frame's usable dynamic range:

```python
row = quality.frame_quality(ext, cal)
print(row["v_sat"], row["v_lim_100sigma"], row["v_lim_5sigma"])
# 8.09  ...  16.40
```

See [depth](depth.md) for the faint end, and note the caveat there — `v_lim_5sigma` runs
optimistic by ~0.7 mag in absolute terms.
