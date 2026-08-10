# Changelog

## Unreleased

## 0.3.0 — 2026-08-10

Aperture sizing is no longer fooled by extended emission, and every photometry entry point
takes a `mask=`. Found while reducing a transit of HD 189733 b, whose field contains M27.

### Fixed

- **Extended emission no longer sizes the aperture.** On a field containing a nebula the
  curve-of-growth cuts — bright, round, isolated — selected the *nebula* rather than stars,
  and the aperture for the whole frame was measured from its growth curve. On a real M27
  600 s stack, of 5469 green detections and 578 bright round sources exactly **one** cleared
  the 40 px isolation cut: M27 itself (`a = 42 px`, `b/a = 0.82`). The resulting green
  aperture was **19.0 px** where the stars wanted 5.3 px, and nothing reported a problem.
  Two compounding causes, both closed:
  - isolation selects *for* extended objects in a crowded field, because only something
    large has no close neighbour, and roundness does not exclude a planetary nebula. New
    `photometry.COG_MAX_SIZE_RATIO` (3.0) rejects candidates whose semi-major axis exceeds
    3× the median of the frame's bright round sources — a star sits at 1, M27 at 21;
  - a COG built from one object returns a *finite* radius, so the documented "falls back to
    `1.2 × FWHM`" never triggered. New `photometry.MIN_COG_STARS` (5) is an explicit floor.

  Only the aperture was affected: masking M27 moved the measured green FWHM by 0.01 px
  (4.41 → 4.40) while moving the 90% aperture from 19.03 to 5.29 px, because `measure_fwhm`
  medians per-source scalars while the COG medians whole curves.
- An empty curve-of-growth sample now returns an explicit all-nan table with
  `meta["n_stars"] = 0` instead of producing the same nans via a numpy "mean of empty slice"
  RuntimeWarning.

### Added

- **`mask=` on the photometry entry points** — `extract_sources`, `forced_photometry`,
  `curve_of_growth`, `measure_fwhm`, `aperture_correction` and `kron.extract_kron` — a
  boolean array, `True` where pixels should be ignored, following SEP's convention. Excludes
  a region from the background estimate, from detection, and from the aperture-sizing
  sample. Two deliberate limits: it is *not* passed to `sep.extract`, because that fragments
  a bright extended source into a ring of detections around the mask boundary whose
  centroids sit outside it (rejecting on the centroid of the un-masked segmentation removes
  the object in one piece); and it never alters a reported flux, so a star whose aperture
  overlaps the mask is still summed over all its pixels.
- **`photometry.sky_mask(shape, wcs, ra, dec, radius_arcsec)`** builds that array from
  circles specified on the sky, which is where the knowledge lives — the pixels move frame
  to frame under dithering and Alt-Az rotation.
- **`Project.mask`** — a `mask(frame, wcs) -> bool array` callable, mirroring the existing
  `provenance` hook, so the batch path can mask per frame; and **`Project.isolation`** to
  relax the curve-of-growth nearest-neighbour requirement. The 40 px default is unreachable
  in a rich field (median nearest-neighbour ~10 px on the M27 fields), which silently forces
  the FWHM fallback and makes `enclosed_lightcurve` inert — check `cog.meta["n_stars"]`.

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
- **ASTAP no longer needs configuring** (`astap`). It is located automatically —
  `$ASTAP_EXE`, then `astap_cli`/`astap` on `PATH`, then a fetched copy, then the stock
  Windows location — so `apt install astap-cli` or any manual install is picked up with no
  `astap_exe=`. `astap.download()` fetches one where there is none: a 1.4 MB binary and a
  100 MB `d05` star database, mirrored under MPL-2.0, unpacked into the cache, made
  executable, de-quarantined on macOS and run once to prove it works.
  `tools/build_astap_mirror.py` assembles the mirror.

### Changed

- `Project.astap_exe` now defaults to `None`, meaning "find it". Passing a path still
  overrides everything. `astrometry.ASTAP_EXE` remains the stock Windows location but is
  now only the last resort rather than the default.

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
