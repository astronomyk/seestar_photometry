# Seestar data formats

Two FITS layouts reach this package, written by different producers, and they disagree
about more than structure. `frames.load_frame` normalises both to one in-memory frame; read
metadata through `frames.frame_metadata`, never straight off a header.

## Layout A — `"cube"`, native Seestar

The telescope writes the debayered, stacked RGB planes as a single 3-D image in the
primary HDU.

```
No.  Name      Type          Dimensions        Format
  0  PRIMARY   PrimaryHDU    (1080, 1920, 3)   int16 (rescales to uint16)
```

FITS axis order is `(NAXIS1=nx, NAXIS2=ny, NAXIS3=3)`, which astropy reads into numpy as
**`(3, ny, nx)`** — channel-first. Verified on real files for both models:

| Model | FITS `(N1, N2, N3)` | numpy shape |
|---|---|---|
| S50 | (1080, 1920, 3) | (3, 1920, 1080) |
| S30pro | (2160, 3840, 3) | (3, 3840, 2160) |
| S30pro, binned | (1080, 1920, 3) | (3, 1920, 1080) |

Note the frames are stored **portrait** (taller than wide). Nothing in the package assumes
otherwise — shape always comes from `frame.shape`.

Stored `uint16`; promoted to `float32` on load so downstream arithmetic cannot overflow.

`load_frame` also transposes a channel-*last* `(ny, nx, 3)` array if it ever sees one.
No Seestar file observed does this, but the failure would be silent and catastrophic —
a colour plane read as image rows — so it is handled rather than assumed away.

## Layout B — `"mef"`, CrowdSky

The CrowdSky platform re-stacks and plate-solves, then writes a multi-extension file:

```
No.  Name        Type          Dimensions      Notes
  0  PRIMARY     PrimaryHDU    ()              empty -- metadata only
  1  RED         ImageHDU      (1080, 1920)
  2  GREEN       ImageHDU      (1080, 1920)
  3  BLUE        ImageHDU      (1080, 1920)
  4  FOOTPRINT   ImageHDU      (1080, 1920)    coverage map, uint8 -- NOT a science plane
  5  STAR-TAB    BinTableHDU   469R x 36C      the server's own SEP + Gaia catalogue
```

Three traps here:

1. **The primary HDU has no `NAXIS1`/`NAXIS2`.** Any code reading image shape from a header
   is broken on these files. Shape must always come from the data.
2. **`FOOTPRINT` is a 2-D image HDU.** Treating it as a science plane would shift the
   colour assignment by one and silently corrupt every colour. Planes are resolved by
   `EXTNAME`, with `FOOTPRINT`/`MASK`/`WEIGHT` explicitly excluded from the positional
   fallback.
3. Frame count varies — some files carry no `FOOTPRINT`, so "5 HDUs" is not a safe test.

`STAR-TAB` is exposed as `frame.star_tab`. It holds the server's SEP detections already
cross-matched to Gaia (`x`, `y`, `flux`, `a`, `b`, `theta`, `flag`, `ra`, `dec`, `gaia_id`,
`gaia_gmag`, `gaia_vmag`, `match_dist_arcsec`). Useful as an independent cross-check; not
used for science, since we re-measure everything.

## Header dialects — and the `EXPTIME` trap

The two producers describe exposure differently, and **`EXPTIME` means different things in
each**:

| Quantity | native Seestar | CrowdSky |
|---|---|---|
| sub-exposure count | `STACKCNT` | `NIMAGES` |
| per-sub exposure | `EXPTIME` (e.g. 10 s) | `EXPTIME / NIMAGES` |
| total on-sky | `TOTALEXP` (e.g. 390 s) | `EXPTIME` (e.g. 410 s) |
| exposure span | `DATE-OBS` only | `OB-START` / `OB-END` |
| wall-clock span | — | `OBSTOTAL` (e.g. 654.9 s) |
| mount mode | `EQMODE` (0 = Alt-Az, 1 = EQ) | absent |

Reading `EXPTIME` naively on a CrowdSky frame gives a **410 s sub-exposure** and, via
`n_exp × exptime`, a ~4.6 hour "total" for a 15-minute stack. The predecessor package did
exactly this: `total_exptime` came back `nan`, so the mid-exposure timestamp defaulted to
`DATE-OBS` and ran **~5 minutes early** — enough to smear and shift a transit.

`frames._exposure` resolves the dialect by *which keywords are present*, not by file
layout, so a future exporter that mixes them still resolves correctly. This is asserted in
`tests/test_frames.py::test_both_dialects_agree_on_physical_quantities`.

### Use the exposure span when it exists

`OB-START`/`OB-END` give the true wall-clock span of the co-add, including the inter-sub
overhead that is otherwise unknowable. `lightcurves.frame_times` prefers it and records
`time_source="span"`. Falling back to `DATE-OBS + total_exptime/2` runs early — for a 390 s
on-sky stack spanning ~650 s of wall clock, by ~2 minutes.

## Model detection

Keyed off **`TELESCOP`**, e.g. `"S50_8a95aa90"` — the model plus the individual unit's
serial. Not `INSTRUME`, which is inconsistent across firmware versions and absent from some
exports. `frames.unit_id` extracts the serial, which matters because each physical Seestar
has its own zero point, PSF and response to conditions.

## The server's own quality block

CrowdSky headers carry a full set of quality metrics from an earlier generation of this
same pipeline — `ZPR/ZPG/ZPB`, `ZPSCT*` (scatter), `ZPCT*` (colour terms), `ZPMAGLO/HI`
(fit window, also 10–14), `ZPNSTAR`, `FWHMR/G/B`, `SKYLVL*`, `SKYRMS*`, `SQMPHOT`,
`BGPOLY0-5`, `BGRESSTD`, `SQMEST`, `BORTEFST`, plus `PLTSOLVD`/`PLTSCALE` from the plate
solve.

`quality.onboard_quality` lifts these into `onboard_*` columns. They are **not** used for
science — everything is re-measured — but comparing them against our own row is a cheap
independent check that a frame reduced sensibly.

## Error handling

The measurement layer does **no** defensive error handling: a broken frame raises. The
batch runner in `pipeline` catches per-frame exceptions, records a status
(`load_error`/`failed`), logs it, and continues. Keeping the two separate is why the science
code stays readable.
