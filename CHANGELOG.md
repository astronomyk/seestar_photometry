# Changelog

## 0.2.0 — 2026-08-05

First release intended for publication.

### Added

- **Raw sub-exposure support.** `load_frame` now recognises the single-plane Bayer layout
  (`BAYERPAT`) and demosaics it on load, so a raw 20 s sub is interchangeable with an
  on-board stack everywhere downstream. The native mosaic is kept on `frame.bayer`.
  New `debayer` module.
- **Local stacking** (`stacking`, needs the `stack` extra): astroalign registration as a
  similarity transform — rotation included, which alt-az frames require — plus a
  footprint-weighted mean. `StackReport` records residual, field rotation, coverage and any
  rejected sub. A local stack writes native-dialect exposure keywords, so `frame_metadata`
  reads it exactly as it reads an on-board stack.
- **Bundled example data** (`examples`): 1000×1000 cutouts of real MW Cam observations —
  three on-board stacks, a CrowdSky multi-extension frame, five raw subs, solved `.wcs`
  sidecars and a trimmed Gaia DR3 table. Every documented example and every real-data test
  runs offline.
- **Documentation** for Read the Docs: install, quickstart, six use-case pages driven by the
  bundled data, an API reference, and the design-decision records.
- 40 real-data tests alongside the existing synthetic suite (162 total).

### Fixed

- `wcs_cache_path` now strips a compression suffix, so `frame.fits.gz` caches its WCS to
  `frame.wcs` rather than the unreachable `frame.fits.wcs`. Gzipped frames would otherwise
  have been re-solved on every run.
- An empty light curve keeps its full column schema, so a downstream call reports "no epochs"
  instead of raising a bare `KeyError: 'dmag'`.
- `pipeline` uses a thread pool when `workers == 1`; a single worker gained nothing from a
  separate process and paid the full spawn cost. A dead process pool now raises a message
  naming the usual Windows cause (a driver without an `if __name__ == "__main__"` guard)
  instead of an opaque `BrokenProcessPool`.
- `forced_photometry` no longer emits a divide warning for a fully off-frame aperture; the
  SNR is `nan`, which is the right answer.
- CrowdSky's `FWHM` header values are in **arcsec** while this package reports **pixels**.
  `quality.onboard_quality` renames them `onboard_fwhm_*_arcsec` so the two are not compared
  by accident — they looked 2.4× apart until the plate scale was applied, after which they
  agree to 2%.

## 0.1.0

Initial extraction of the photometry engine from the MW Cam analysis repo, where it had
existed as three diverging copies. See `docs/migration-from-mwcam.md`.

Verified against the frozen results it replaced: zero point bit-identical on 20 of 22 frames
(all 15 photometric-grade frames exactly 0.000 mmag apart), and forced photometry
bit-identical across a 52-frame CrowdSky set.
