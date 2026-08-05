# Quickstart

Everything on this page runs offline against the bundled example frames.

## Load a frame

```python
from seestar_photometry import examples, frames

frame = examples.stack()          # or: frames.load_frame("my_stack.fit")
print(frame.layout, frame.shape, frame.model)
# cube (1000, 1000) S50
```

Whatever the file layout — on-board stack, raw Bayer sub, or a CrowdSky multi-extension
frame — you get the same thing: a `(3, ny, nx)` float32 cube in R, G, B order.
**Green is index 1 and is the science band** (≈ Johnson V through the IRCUT filter).

```python
frame.g          # the green plane, (1000, 1000)
frame.shape      # (ny, nx) -- always from the data, never a header
```

## Read the metadata

Never read exposure keywords off the header directly. `EXPTIME` means the *per-sub*
exposure in a native file and the *total* on-sky time in a CrowdSky one; taking it at face
value gives a mid-exposure timestamp five minutes wrong. `frame_metadata` resolves the
dialect:

```python
meta = frames.frame_metadata(frame)
print(meta["n_exp"], meta["exptime"], meta["total_exptime"], meta["unit"])
# 38 20.0 760.0 98ac1c17
```

38 sub-exposures of 20 s, 760 s on sky, from the Seestar with serial `98ac1c17`.
`airmass` is computed from the pointing, site and time — Seestar headers carry none.

## Detect and measure

```python
from seestar_photometry import photometry

ext = photometry.extract_sources(frame)
print(len(ext.band("G")), "stars above SNR 5 in green")
# 362 stars above SNR 5 in green
```

The aperture is sized **per frame and per band** from each plane's curve of growth, to the
radius enclosing 90% of the flux:

```python
print(ext.fwhm)       # [4.35 3.80 4.65]  px, R G B
print(ext.aperture)   # [4.67 4.37 6.35]  px, R G B
```

Red is 1.15× and blue 1.22× the green FWHM. That is the Seestar's chromatic aberration, not
noise, and it is why a single shared radius would bias every colour — see
[photometry-design](photometry-design.md).

## Calibrate onto Gaia V

```python
from seestar_photometry import calibration

ext.match_gaia(examples.gaia(), wcs=examples.wcs())
cal = calibration.fit_zeropoint(ext.sources, band="G")
print(cal)
```

```
zeropoint   = 23.644 ± 0.003 mag
colour_term = +0.036   at colour0 = 1.137 (B-R, JKC)
rms         = 0.016 mag over 59 stars
```

The fit is `V = m_inst + ZP + k·(B−R − colour0)`, restricted to V ∈ [10, 14], with
catalogue variables and sigma-clipped outliers excluded. An rms of 0.016 mag is a good
frame; below 0.06 counts as photometric-grade.

```{admonition} chi2_red will look alarming
:class: note
This fit reports `chi2_red = 214`. That is normal. It measures scatter against the
*photon-noise* prediction, and Seestar frames sit on a ~0.03 mag systematic floor. Use it to
spot frames far worse than their peers, not as an absolute goodness test.
```

## Where did that WCS come from?

The bundled frames ship with solved `.wcs` sidecars. For your own frames you must solve
once — the on-board WCS is off by ~1 arcminute, which drops the match rate to a few percent:

```python
from seestar_photometry import astrometry

wcs = astrometry.solve(frame, solver="astap")   # cached beside the FITS; ~1 s
wcs = astrometry.load_wcs(frame)                # reads the cache on later runs
```

Check it worked by the only measure that matters — how well the sources land on the
catalogue:

```python
print(astrometry.match_quality(ext.band("G")))
# {'median_arcsec': 0.51, 'p90_arcsec': 2.50, 'matched_frac': 0.88, 'n': 362}
```

A 0.51″ median with 88% matched is a good solve. A median approaching the 2″ tolerance means
the solve is wrong even though it "succeeded".

## Per-frame quality

```python
from seestar_photometry import quality

row = quality.frame_quality(ext, cal)
print(row["v_lim_5sigma"], row["sky_sb"], row["bortle"])
# 17.15 18.72 7
```

V ≈ 17.15 at 5σ in 760 s, under a Bortle 7 sky. That is one row of the frame table
described in [frame-table](frame-table.md); `pipeline` builds them in bulk.

## See it, don't trust it

```python
from seestar_photometry import report

report.frame_report(frame, ext, cal, "diagnostics/")
```

Thirteen figures: curve of growth with the chosen aperture marked, per-band FWHM ratios, the
zero-point relation and colour-term fit, residuals against magnitude, SNR and position, the
background triptych, detections with apertures drawn on, and the cross-match separation
histogram. [diagnostics](diagnostics.md) says what each one answers and what wrong looks
like.

## Next

- [Stacking raw subs](usecases/stacking.md) — from 20 s Bayer frames to a stack
- [Zero point](usecases/zeropoint.md) — calibration in detail
- [Saturation limit](usecases/saturation.md) — the bright end
- [Depth vs exposure](usecases/depth.md) — how deep, and how deep you *could* go
- [A light curve](usecases/lightcurve.md) — ensemble differential photometry
- [Batch processing](usecases/batch.md) — hundreds of frames, resumably
