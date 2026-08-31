# Handoff: all-sky synthetic R_S50 / G_S50 / B_S50 magnitudes for the gaiadb

*From the Claude instance in `photometry-adventures-2`, 2026-08-31. That repo
derived and validated the Seestar S50's native photometric system; this job
scales the synthesis from one field (13k stars) to the full catalogue this
package hosts. Planned execution host: **Zeus** (32 cores / 256 GB RAM / 6 TB
SSD / fibre), since the local workstation is commit-limited.*

## Goal

Add synthetic S50-system magnitudes — `r_s50_mag`, `g_s50_mag`, `b_s50_mag`
(+ `_err`, AB scale) — to every row of the `gaia-seestar-dr3-v1` dataset
(the HEALPix-5-partitioned Parquet catalogue built by `tools/build_gaia_allsky.py`,
GSPC sources cut at V < 17.5). Per `gaiadb.py`'s own convention this is a schema
change: bump `DATA_VERSION` to `v2`.

Why: calibrating the Bayer planes against these references instead of
`v_jkc_mag`/`b_jkc_mag`/`r_jkc_mag` removes the ZP colour term entirely.
Measured on 57 fifteen-minute 20 Cep stacks: G +0.054 → −0.002, R −0.132 →
−0.002, B +0.548 → ≈0 mag/mag; B-plane calibrator rms 347 → 104 mmag. Full
story: `photometry-adventures-2/docs/handoff-bayer-synthetic-photometry.md`
and the write-up artifact
https://claude.ai/code/artifact/b6dbd950-1df7-4140-9bea-24634374277a

## Inputs (must travel to Zeus)

1. **The passband curves** — copy, do not rebuild:
   `D:\CrowdSky\crowdsky-workspace\Repos\photometry-adventures-2\paper\data\bayer_passbands.ecsv`
   Columns `wavelength_nm, R_S50, G_S50, B_S50` (relative response, 390–1000 nm,
   2 nm). Model: IMX462 datasheet curves × logistic UV/IR-cut (half power
   390 / 693 nm, 10–90% width 15 nm) with the B dye's red tail scaled
   **×1.15** above 540 nm. **The adopted tail is 1.15** (the three-frame scan
   said 1.18, the full 57-stack survey preferred ~1.12; 1.15 is the adopted
   compromise) — verify the ECSV meta says 1.15 before shipping; if it says
   anything else, rebuild via `paper/scripts/build_bayer_passbands.py` in that
   repo first.
2. **The reference implementation**:
   `photometry-adventures-2/paper/scripts/gaia_bayer_synth.py`. Two functions
   matter: `synth_mags()` (photon-weighted AB integration — copy it verbatim,
   including the `resp * wl` photon weighting and `np.trapezoid`) and
   `load_spectra()`/`bulk_files_for()` (the bulk fetch, below).
3. **The gaiadb source_id sets** — either the parquet dataset itself or just
   the per-part `source_id` columns; the job only needs ids in, mags out.

## Data source: ESA bulk CDN (NOT the archive query services)

The Gaia TAP/datalink services hung for 15+ minutes on 100-source requests on
two separate evenings. The bulk product is served from a fast CDN and is the
only sane route at this scale:

- Product: `xp_sampled_mean_spectrum` — calibrated spectra, `flux` and
  `flux_error` as float32[343] on a fixed **336–1020 nm, 2 nm** grid, units
  W m⁻² nm⁻¹. Gzipped ECSV, one file per HEALPix-8 range, named
  `XpSampledMeanSpectrum_<hp8lo>-<hp8hi>.csv.gz`.
- Download base: `https://cdn.gea.esac.esa.int/Gaia/gdr3/Spectroscopy/xp_sampled_mean_spectrum/`
- Enumeration: the S3-style listing at `https://gaia.eu-1.cdn77-storage.com/`
  with `?prefix=Gaia/gdr3/Spectroscopy/xp_sampled_mean_spectrum/XpSampledMeanSpectrum_`
  and `&marker=<last key>` pagination (`<Key>`/`<Size>`/`<IsTruncated>` XML).
- **Enumerate first, then commit**: list every file with sizes and report
  count × total GB before downloading. The six files covering the 20 Cep field
  were 30–41 MB / 9k–12k rows each; scale up from the real listing, not from
  guesses. Zeus's 6 TB absorbs any plausible total; keep the files after the
  job — they are reusable for any future field or re-synthesis.

## The id ↔ file mapping (this is the whole trick)

A Gaia DR3 `source_id` encodes its nested-HEALPix position in the high bits:

    hp5 = source_id >> 49      # the gaiadb parquet partition (gaiadb.hpx_of)
    hp8 = source_id >> 43      # the bulk-file partition key

Each hp5 pixel is exactly 64 consecutive hp8 pixels (`hp5*64 .. hp5*64+63`),
so every parquet part maps to a contiguous hp8 range and therefore to a short
contiguous run of bulk files. Recommended shape, **driven by bulk file** so
each is read exactly once:

    for each bulk file (parallel over files):
        rows = parse(file)                       # source_id, flux, flux_error
        for each hp5 part its hp8 range touches:
            keep = rows whose source_id is in that part's id set
            integrate keep through the three passbands -> mags
            append (source_id, 6 columns) to that part's staging output
    then per part: join staging onto the parquet part, write v2 part

## Numerics

    f_AB(λ)  = 3631e-26 W m⁻² Hz⁻¹ × c / λ²          (c = 2.99792458e17 nm/s)
    weight   = S(λ) × λ                                (photon counting)
    mag      = −2.5 log10( ∫ f_λ w dλ / ∫ f_AB w dλ )
    mag_err  = 1.0857 × sqrt(Σ (σ_λ w Δλ)²) / ∫ f_λ w dλ

Interpolate each passband onto the file's 336–1020 grid with zeros outside
390–1000. The error ignores XP sample correlations (mild underestimate —
documented, accepted). Round mags to 4 decimals.

## Performance notes

- Local timing: ~5–8 s per 35 MB file for astropy-ECSV parse + filter +
  integrate, and the **parse dominates**. At bulk scale consider skipping
  `Table.read`: strip the ECSV header, feed the CSV body to pyarrow/pandas
  with array-string converters — but validate against astropy output first.
- Parallelize over files, ~24–28 workers on Zeus; the integration is trivial
  vectorized numpy. Set `OPENBLAS_NUM_THREADS=1` (etc.) per worker.
- Not every gaiadb row has an XP spectrum (20 Cep: 12,977 of 13,086). Missing
  → NaN in the six new columns; never drop rows. Report global coverage.

## Validation (required before adopting)

1. Extract the 13,086 20 Cep source_ids from
   `photometry-adventures-2/paper/data/gaia_20Cep.ecsv` and compare the job's
   mags against `paper/data/gaia_20Cep_bayer.ecsv` — but note that file is
   whatever tail scale that repo last built; regenerate it at 1.15 there if in
   doubt. Agreement should be ≤ 0.0001 mag (identical inputs and formula).
2. Spot-check a sparse high-latitude part and a Galactic-plane part for
   coverage and sane colour distributions (G_S50−R_S50 roughly 0.1–0.5 for
   most stars).

## Naming

Bands are named per unit model: `*_s50`. The S30, S30 Pro and S50 Pro use
different sensors (Sony IMX585, OmniVision OS08B10) and will get their own
fitted curves (`*_s30`, ...) once measured — leave room in the schema
convention for those, don't generalize the column names.
