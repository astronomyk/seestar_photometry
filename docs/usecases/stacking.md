# Stacking raw sub-exposures

The Seestar stacks on board, but only into its own bins, and only after every sub has made
the round trip to the scope and back. Stacking locally from the raws lets you choose the bin
boundaries, keep every sub, and stack the same photons more than one way.

Needs the `stack` extra: `pip install "seestar-photometry[stack]"`.

## A raw sub is not three planes

An on-board stack arrives debayered. A raw sub does not — it is a single 2-D mosaic tagged
by `BAYERPAT`:

```python
from seestar_photometry import examples, debayer

raw = examples.raw_subs(1)[0]
print(raw.layout, raw.shape, raw.bayer.shape, debayer.pattern_of(raw.header))
# bayer (1000, 1000) (1000, 1000) GRBG
```

`load_frame` demosaics it for you, so `raw.data` is the usual `(3, ny, nx)` cube and a raw
sub is interchangeable with a stack everywhere downstream. The undemosaiced mosaic stays on
`.bayer`, because interpolated pixels are correlated and both registration and any honest
per-pixel noise estimate want the native samples.

The channel balance of a raw sub is worth seeing once:

```python
print(debayer.channel_medians(raw.bayer))
# {'G': 965.0, 'R': 986.0, 'B': 1103.0}
```

Blue reads ~140 ADU above green. The on-board stacker balances these to ~963/965/964; a
local stack does not. It makes no difference to the photometry, which lives entirely in
green.

```{admonition} The mosaic phase is not a guess
:class: note
For `GRBG` as the array is indexed, green sits on `(0,0)` and `(1,1)` of each 2×2 quad, red
on `(0,1)`, blue on `(1,0)`, with no row flip. The check: the two green sub-lattices agree
to 4 ADU with matching sigma, while `(0,1)` reads 986 and `(1,0)` 1103. Get the phase wrong
and you silently swap colour planes.
```

## Co-add

```python
from seestar_photometry import stacking

frame, report = stacking.stack_frame(examples.raw_sub_paths())
print(report)
# 5/5 subs stacked (100 s on sky), residual 0.13 px, field rotation 0.31 deg, coverage 100.0%
```

Read the report before trusting the result. Three numbers matter:

`residual 0.13 px`
: Median registration residual. Sub-pixel is what you want. A sub whose fit is worse than
  1 px, or which matches fewer than 6 sources, is **rejected** rather than averaged in — a
  bad transform smears every star in the stack.

`field rotation 0.31 deg`
: Across 100 s. These are alt-azimuth frames, so the field genuinely rotates; over a
  15-minute bin it reaches several degrees. Registration fits a *similarity* transform
  (shift, rotation, scale) for exactly this reason. A rotation span of ~0 would mean the fit
  collapsed to a pure shift, which is the failure that smears the frame corners.

`coverage 100.0%`
: Fraction of output pixels any sub reached. Rotation leaves corners uncovered; those pixels
  are zero and excluded by the weight map, rather than averaged against zeros and dragging
  the sky low.

`report.reasons` lists any rejected sub and why.

## Did it help?

```python
from seestar_photometry import photometry

one = photometry.extract_sources(examples.raw_subs(1)[0])
many = photometry.extract_sources(frame)

print(f"green sky RMS: {one.rms[1]:.1f} -> {many.rms[1]:.1f} ADU")
print(f"stars above SNR 5: {len(one.band('G'))} -> {len(many.band('G'))}")
# green sky RMS: 41.9 -> 17.6 ADU
# stars above SNR 5: 264 -> 489
```

The noise falls by 2.38× on five subs, against √5 = 2.24 — slightly better than the ideal,
because the weighted mean plus demosaic correlation buys a little extra smoothing. The
detected-source count nearly doubles.

## The stack reads like an on-board one

```python
from seestar_photometry import frames

print(frames.frame_metadata(frame))
# n_exp=5, exptime=20.0, total_exptime=100.0
```

`STACKCNT` and `TOTALEXP` are written in the native Seestar dialect deliberately, so
`frame_metadata` reads a local stack exactly as it reads an on-board one. Without that the
per-sub 20 s would be mistaken for the total and every mid-exposure timestamp would be wrong.

The reference sub's WCS is **dropped**: it no longer describes the co-add, and leaving it in
place would look solved while reintroducing the ~1 arcminute error the package re-solves to
avoid. Solve the stack like any other frame.

```python
stacking.write_stack("my_stack.fits", frame.data, frame.header, report)
```

Writes a single `(3, ny, nx)` primary HDU — the on-board layout, indistinguishable to
anything reading it.

## What this deliberately does not do

No outlier rejection, no per-frame quality weighting, no channel balancing, no gradient
removal. The omissions are the point: they make the local stack a fair test of whether the
on-board stacker is doing anything a plain pipeline cannot. On the MW Cam data it was not —
there was no photometric advantage either way, which means the raws are worth archiving and
the round trip is optional.

For thousands of subs, `coadd` streams one frame at a time into a single accumulator, so
memory stays flat.
