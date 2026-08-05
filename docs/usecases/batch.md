# Batch processing a whole dataset

A night is hundreds of frames and a season is thousands. Three stages, all resumable.

## Configure once

There are no config files and no CLI. A `Project` is the entire configuration surface, so
retargeting means constructing a different object rather than editing code:

```python
from seestar_photometry import LocalTree, Project, Target, pipeline

proj = Project(
    target   = Target("MW Cam", ra=186.6821, dec=81.474, source_id=1719678555695057280),
    source   = LocalTree(roots=[r"D:\data\MW Cam\stacks"]),
    work_dir = r"D:\work\mwcam",
    solver   = "astap",
)
```

Set `source_id` when you know it. Without it the target is taken to be the catalogue source
nearest the pointing, which is right for an isolated centred target and wrong in a crowded
field.

## Three stages, in order

```python
pipeline.solve_all(proj)                                  # .wcs sidecars
frames = pipeline.build_frame_table(proj, diagnostics=3)   # frames.ecsv
stars, meas = pipeline.build_measurements(proj)            # stars.ecsv, measurements.ecsv
```

The order is required: stages 2 and 3 read the cached WCS and never solve. That split is
deliberate — solving is network- or subprocess-bound and its cache is per-frame and reusable,
while characterisation and measurement are pure CPU once the WCS and catalogue are cached, so
they parallelise across cores and can be re-run freely as choices change.

```python
pipeline.build_frame_table(proj, workers=8)
```

## Resumability is a normal operation

Interrupt any stage and re-run it. Frames already in the output table, or already carrying a
sidecar, are skipped; rows are checkpointed to disk every 25 frames. Killing a run that has
hours left and restarting later is expected usage, not a recovery procedure.

A single bad frame never sinks a batch. Each is caught and recorded with a status:

| Status | Meaning |
|---|---|
| `ok` | measured |
| `cached` | already done, skipped |
| `no_wcs` | no sidecar — run `solve_all` first |
| `load_error` | unreadable FITS (a truncated download) |
| `failed` | raised during measurement; logged to `work_dir/errors.log` |

## What you get

Everything lands in `work_dir` — never in the data tree, never in the FITS headers. The one
exception is the per-frame `.wcs` sidecar, which lives beside its frame because it is
expensive to recompute and useful to every project touching that frame.

```python
frames, stars, meas = pipeline.load_tables(proj)
```

On the bundled three-frame example:

```
frames.ecsv    3 rows, 57 columns
  stack_saturated   ZP=23.568  rms=0.051   280 s
  stack_c17_15min   ZP=23.644  rms=0.016   760 s
  stack_c17_30min   ZP=23.656  rms=0.017  1460 s
```

`frames.ecsv` is the table you query to decide what to trust — see
[frame-table](../frame-table.md).

```python
from seestar_photometry import quality
good = quality.photometric(frames)      # rms < 0.06
```

## Project-specific columns

Datasets carry bookkeeping the package should not know about — a dataset name, a stacking
manifest, binning. Supply it through a hook rather than forking the row builder:

```python
def provenance(frame):
    return {"dataset": frame.path.parent.parent.name, "observer": "JDB"}

proj = Project(..., provenance=provenance)
```

It must be a module-level function: the batch runner ships it to worker processes.

## Curation is yours, deliberately

Which frames of *your* dataset to use is a judgement the package will not make for you —
baking a rule in would silently apply it to unrelated data. `LocalTree` takes a predicate:

```python
from seestar_photometry.frames import object_name

def mine(path):
    obj = object_name(path)              # reads one header, cheaply
    return obj == "" or "MW" in obj      # a blocklist, not an allowlist

proj = Project(source=LocalTree(roots=[...], curate=mine), ...)
```

A blocklist is usually right: object labels are inconsistent across contributors and some
units write no `OBJECT` at all, so keeping anything plausible and dropping only what
positively names another target loses fewer good frames.

## Two aperture settings, on purpose

```python
proj.enclosed_characterise   # 0.90 -- frame table, depth limits
proj.enclosed_lightcurve     # 0.95 -- forced photometry for light curves
```

They optimise different things and should not be unified. 0.90 sits near the SNR sweet spot
and is what published depth numbers use; 0.95 minimises the *differential* scatter floor,
because a slightly larger aperture damps the position-dependent aperture loss from alt-az
field rotation. See [photometry-design](../photometry-design.md).

## A Windows footgun

```python
if __name__ == "__main__":
    main()
```

Wrap your driver. On Windows, `ProcessPoolExecutor` re-imports the main module in every
worker, so a script that runs a stage at import time spawns endlessly and the pool dies. The
package detects this and raises a message saying so, but the guard is the fix. Or pass
`workers=1`, which uses threads and skips process spawning entirely.
