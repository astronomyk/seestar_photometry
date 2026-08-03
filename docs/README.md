# Documentation

These are **decision records**: *why* the code is the way it is, with the measurements
behind each choice. They are not API reference — the docstrings are.

Read them when you are about to change a numerical choice, or when a result looks wrong and
you need to know what was already ruled out.

| Doc | Topic |
|---|---|
| [photometry-design.md](photometry-design.md) | **Start here.** Fixed circular aperture over Kron; aperture sized per frame *and per band* by enclosed flux; why two different fractions; the aperture correction |
| [data-format.md](data-format.md) | The two FITS layouts, the two header dialects, and the `EXPTIME` trap |
| [astrometry-and-gaia.md](astrometry-and-gaia.md) | Why the on-board WCS is unusable; ASTAP vs astrometry.net; the cached Gaia mosaic; the zero-point fit |
| [frame-table.md](frame-table.md) | The `frames.ecsv` schema, how to read `rms` and `chi2_red`, why `v_lim_5sigma` is optimistic, the condition-corrected depth model |
| [light-curves.md](light-curves.md) | The two-table model; forced photometry; timing; choosing comparisons; the ensemble differential and why it is robust |
| [diagnostics.md](diagnostics.md) | Every figure, what it answers, and what wrong looks like in it |
| [architecture.md](architecture.md) | Module layers, the conventions and their reasons, the batch runner, the frame-source seam |
| [migration-from-mwcam.md](migration-from-mwcam.md) | What changed from the predecessor package, and how to port a driver |

`../CLAUDE.md` is the short version of the conventions, for when you are editing code.

## The three findings that most shape this package

1. **The Seestar PSF is chromatic** — R and B are measurably broader than green, and the
   broadening tracks auto-focus. Hence a per-band aperture sized to a common *fraction*
   rather than a common radius. → [photometry-design.md](photometry-design.md)
2. **The on-board WCS is off by ~1 arcmin**, which drops the catalogue match rate to a few
   percent. Every frame is re-solved and cached. → [astrometry-and-gaia.md](astrometry-and-gaia.md)
3. **Each comparison must be referenced to its own catalogue magnitude**, not to an ensemble
   mean flux. That is what makes a changing ensemble harmless and makes the per-comparison
   scatter a real error estimate. → [light-curves.md](light-curves.md)
