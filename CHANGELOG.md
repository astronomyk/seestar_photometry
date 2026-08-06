# Changelog

## Unreleased

### Added

- **A catalogue-anchored plate solver** (`solver="local"`). No external binary, no network
  and no index files: it pairs SEP detections against the reference catalogue the pipeline
  already needs. Not a blind solver, because it does not need to be — the header pointing
  says roughly where the telescope looked, so the job is to pin the field down rather than
  find it. It seeds from the header WCS where there is one and otherwise votes on
  detection-minus-catalogue offsets over a rotation grid, then fits least-squares over every
  matched pair. Validated against ASTAP on the MW Cam datasets: the two agree to a median
  0.47 arcsec on S50 raw subs, 1.12 arcsec on S30pro subs, and 0.30/2.05 arcsec on locally
  built 600 s stacks — all well under a pixel. `solver="astap"` remains the default: it is
  several times faster (0.2–0.4 s against 1–7 s) and needs no catalogue.
- **An offline Gaia catalogue** (`gaiadb`, needs the `catalog` extra): an opt-in download
  that removes the TAP query, the least reliable step in the pipeline. HEALPix level-5
  partitioned Parquet, so a region can be fetched on its own — three fields cost tens of MB
  rather than the whole multi-GB set. `Project(catalogue_backend="auto")` uses it when it
  covers the field and falls back to TAP when it does not, writing the same ECSV either way,
  so nothing downstream changes.
- **Proper motions and epoch propagation.** `Project(epoch=2026.4)` moves catalogue
  positions from Gaia's J2016.0 to the observing season. A decade of proper motion moves a
  100 mas/yr star by 1 arcsec — half the default match tolerance. The TAP path never fetched
  proper motions and so could not do this at all.
- The offline catalogue also carries columns nothing reads yet but that are cheap and hard
  to add later: `c_star` and `ipd_frac_multi_peak` (blend indicators), `ruwe` and
  `non_single_star` (unresolved binaries), `teff_gspphot`, and `v_jkc_flag`.
- `tools/build_gaia_catalogue.py` builds a conforming dataset for one region from TAP,
  which is enough to use the offline path on a field before any full-sky build exists.

### Fixed

- `docs/astrometry-and-gaia.md` claimed every Seestar frame arrives with a header WCS.
  Native stacks carry the pointing and the optics keywords but no `CRVAL` at all; only
  CrowdSky frames have a header solution.

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
