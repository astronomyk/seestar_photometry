# The example data

Real Seestar observations back every example in these docs and every real-data test. They
are fetched on demand rather than bundled, and cached, so the wheel stays tiny:

```python
from seestar_photometry import examples

examples.download()          # ~18 MB, once; every accessor calls it implicitly
print(examples.available())
print(examples.data_dir())   # where it went
```

See [install](install.md#the-example-data) for the cache location and the offline route.

| Accessor | File | What it is |
|---|---|---|
| `examples.stack()` | `stack_c17_15min` | 760 s c17 stack, rms 0.016 mag, 59 calibrators. The workhorse. |
| `examples.stack_deep()` | `stack_c17_30min` | 1460 s stack of the same field, two hours later |
| `examples.stack_saturated()` | `stack_saturated` | 280 s stack whose brightest star clips (V ≈ 8.1) |
| `examples.crowdsky()` | `crowdsky_mef` | CrowdSky multi-extension frame, plate-solved server-side |
| `examples.raw_subs()` | `raw_sub_1..5` | Five consecutive 20 s raw Bayer subs |
| `examples.gaia()` | `gaia_mwcam.ecsv` | 830 Gaia DR3 rows with synthetic Johnson V |
| `examples.wcs(name)` | `*.wcs` | The solved WCS of a bundled frame |

Convenience helpers:

```python
examples.target()             # a Target for MW Cam
examples.project("work/")     # a ready-to-run Project, Gaia cache seeded
examples.raw_sub_paths()      # paths, in time order -- what stacking wants
examples.path("gaia_mwcam")   # any example file by stem
```

## What they are, exactly

Every frame is a **1000×1000 cutout** of a real observation, gzipped — about 18 MB in total.
Real pixels, real headers, real stars. Full frames are ~15 MB *each*.

The archive is a GitHub release asset (`example-data-v1`), verified against a pinned SHA-256
on download and unpacked through a staging directory, so an interrupted fetch cannot leave a
half-populated cache that then looks complete.

The WCS stays valid because the cutout origin is subtracted from `CRPIX`; nothing else in the
solution changes, since plate scale, rotation and SIP distortion are properties of the optics
rather than of which sub-region was kept. Each frame's `SPCUTOUT` and `SPORIGIN` header cards
record the region and the source file.

All of them show the same field — **MW Cam**, RA 186.6821, Dec +81.474 — so one Gaia table
serves everything. Being at Dec +81 is incidentally useful: it is where naive RA arithmetic
breaks, so the examples exercise the awkward case.

The three stacks are chosen from a per-frame quality table rather than at random:

- the two c17 frames are the same telescope on the same night at two integrations, which is
  what makes a zero-point cross-check and a depth comparison meaningful;
- the saturated frame's cutout is deliberately **centred on its clipping star**, because a
  centred cutout dropped that star and `saturation_mag` then returned `nan` — which is
  precisely what the example needs to demonstrate.

## Limits

The dataset is chosen to exercise code paths honestly, not to reproduce a science result:

- **Three epochs is not a light curve.** The mechanics work, but no signal can be recovered
  and the per-comparison diagnostic needs tens of epochs to mean anything. The full MW Cam
  numbers are quoted in [the light-curve use case](usecases/lightcurve.md).
- **A blind plate solve cannot be demonstrated** on a 0.35° cutout. The solved sidecars ship
  instead; solving is shown against your own full frames.
- **Only S50 frames.** The S30pro's 3840×2160 frames are four times the pixels and would
  dominate the wheel.

## Regenerating

The full datasets are not in the repo, so this only works on a machine that has them:

```bash
uv run python tools/build_example_data.py --size 1000
```

It writes the files to `example_data/`, packs them into `example-data-v1.tar.gz` and prints
the SHA-256 to paste into `examples.DATA_SHA256`, then tells you the `gh release create`
command. The script documents which source frame each example comes from and why.

Bumping `DATA_VERSION` in both the tool and `examples.py` is what invalidates every cache —
no one has to clear anything by hand.
