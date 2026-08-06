# Astrometry and the reference catalogue

## The on-board WCS is unusable

What a frame arrives with depends on who wrote it, and neither case is good enough:

- A **native Seestar stack** carries the pointing (`RA`/`DEC`) and the optics keywords, but
  no WCS at all — verified on every bundled example frame, none of which has a `CRVAL`.
- Where a header WCS *is* present it is **not** photometric-grade: positions are wrong by
  tens of pixels (~1 arcmin), which drops the catalogue match rate to a few percent and
  makes forced photometry meaningless — the aperture lands on empty sky.

The exception is CrowdSky, which plate-solves server-side and says so with `PLTSOLVD`.

So every frame is re-solved, and the solution is cached as a `.wcs` sidecar (a header-only
FITS) **next to the FITS file**. That location is deliberate and is the one exception to
"derived products go in `work_dir`": a solve is expensive, is a property of the frame rather
than of any analysis, and is reused by every project that touches that frame.

A good solve gives a **median cross-match separation of ~0.8 arcsec**. The
`*_match_sep` diagnostic panel is how you check: a distribution filling the match tolerance,
or pressed against it, means the solve is wrong even though it "succeeded".

## Three solvers, one sidecar

`astrometry.solve(frame, solver=...)`:

**`"local"`** — anchored on the reference catalogue rather than on an index of quads, and
**not a blind solver, because it does not need to be.** The pointing tells you roughly where
the telescope was looking, so the field is already known; what is left is to pin it down.
That makes it the only solver with no external dependency of any kind — no binary, no
network, no index files.

Two routes to an approximate solution, tried in order:

1. **The header WCS**, when there is one — which in practice means CrowdSky, since native
   Seestar frames carry none. It is used as a *starting point*, never trusted; that is what
   separates it from `"lift"`.
2. **Voting**, otherwise. The orientation is unknown and the pointing is off by tens of
   arcmin, but the plate scale is known to about 1%. So for each angle on a 2° grid, the
   catalogue's predicted pixel positions are rotated and *every* detection-minus-prediction
   offset is histogrammed. Shared stars pile into one bin; everything else spreads flat. The
   sharpest peak gives the orientation and the shift, and averaging the pairs in the winning
   bin gives the shift to well under a pixel. A second, one-dimensional scan then recovers
   the plate scale.

Either way the approximate solution is only a seed: the returned WCS is a least-squares fit
(`astropy.wcs.utils.fit_wcs_from_points`) over every matched pair, iterated with a tightening
tolerance, and rejected outright if too few pairs survive.

### Why not an asterism matcher

The obvious choice was astroalign's triangle matcher, which `stacking` already uses between
raw subs. It works on the bundled cutouts and **fails on every real full frame**, and the
reason is worth recording so it is not tried again.

Triangle invariants are scale-free, so the matcher throws away the one thing we know for
certain and has to rediscover it. It can also only match the brightest N of each list — and
the catalogue has to cover the frame's diagonal *plus* the pointing error, which for an S50
sub is 6.2 deg² against a 0.92 deg² frame. Only ~15% of the brightest catalogue sources are
on the frame at all, so the two brightness rankings simply disagree. No amount of tuning
`max_control_points` fixes that; it only makes the failure slower. Voting has no such
problem, because contamination adds a flat background rather than competing with the signal.

### The things that bite

> **The Seestar image is mirrored.** Its solved CD matrix has a *positive* determinant:
> right ascension increases with *x*, the opposite of the usual north-up/east-left
> convention. Measured on both an S50 and an S30pro. A rotation cannot express a reflection,
> so seeding with the wrong parity does not merely slow the search down — it makes it
> impossible. `astrometry.SEED_PARITY` tries the mirrored orientation first and the
> conventional one second, and stops as soon as one wins decisively.

**The pointing is off by far more than an arcmin.** Measured against ASTAP on raw subs: a
remarkably constant 9.0 arcmin on every S50 frame tried, and 9.3 to 40.1 arcmin on the
S30pro. The S50's constancy suggests a fixed offset between the reported pointing and the
sensor centre rather than mount error. `POINTING_SLACK_ARCMIN` sets how far the search
reaches.

**Field rotation is large.** A bundled 15-minute c17 stack needed **33 degrees**. Nothing
shift-only would work.

**The nominal plate scale is not good enough to pair with.** `206.265 × XPIXSZ / FOCALLEN`
ran 0.8% high on the S50 and 1.5% high on the S30pro against the solved value. On an
S30pro's 2200-pixel half-diagonal that 1.5% is 33 pixels — further than the typical distance
to the *wrong* neighbour, so pairing at the field edge would silently pick the wrong star.
Hence the separate scale scan, and hence the rotation vote using only the core of the frame,
where the error is small.

**Solving detects at 5σ, not the photometry's 2σ.** `Project.thresh` is kept low on purpose
so that the SNR cut and not the detection defines the sample. A solver wants the opposite:
faint detections are below the reference catalogue's limit, so they cannot pair with
anything and only offer more ways to go wrong. `detect_for_solve` also skips
`extract_sources` entirely — apertures, curves of growth and forced photometry have no
bearing on where the stars are, and on a deep S30pro co-add they cost **186 s** against 8 s
for the detection alone.

### Measured performance

Against ASTAP on the MW Cam datasets, worst-case disagreement over the frame and median
per-frame wall clock:

| | frames | solved | vs ASTAP (median) | local | ASTAP |
|---|---|---|---|---|---|
| S50 raw subs (20 s) | 8 | 8 | 0.47″ | 1.0 s | 0.18 s |
| S30pro raw subs (30 s) | 8 | 6 | 1.12″ | 2.9 s | 0.36 s |
| S50 stack (600 s, local) | 1 | 1 | 0.30″ | 2.8 s | 0.16 s |
| S30pro stack (600 s, local) | 1 | 1 | 2.05″ | 7.2 s | 0.37 s |
| bundled example stacks | 4 | 4 | 0.09″ | 0.5 s | — |

**ASTAP is several times faster and does not need a catalogue.** It stays the default. The
local solver earns its place by needing nothing installed, and by being accurate enough that
the choice is about convenience rather than quality — on the bundled stacks its cross-match
medians (0.41–0.49″) beat ASTAP's 0.51″ on the same frames.

The two S30pro subs that did not solve failed *loudly*, which is the point: a returned
solution has to be trustworthy. `MIN_PAIR_FRACTION` is what enforces that — a wrong solve
pairs almost nothing proportionally (measured at 2% on a sub that came out 65″ off) against
18–55% for every correct one, including a cloud-affected frame where most detections are
spurious. An absolute floor is not enough on a rich field, where a few thousand detections
against a few thousand catalogue sources throw up coincidences whatever the WCS says.

**`"astap"` (default)** — a local, fully offline blind plate solver. Nothing needs
configuring: `astap.executable` looks for `astap_cli` or `astap` on `PATH` before falling
back to a copy `astap.download()` fetched, and only then to the stock Windows location. The
old default was a hardcoded `C:\Program Files` path, which meant every non-Windows user had
to set `astap_exe=` by hand even after `apt install astap-cli`.

There is no PyPI package that bundles ASTAP — `astapy` wraps the CLI but still expects a
manual install, which is the problem rather than the solution — so the fetch is a mirror of
the upstream binaries. ASTAP is MPL-2.0, so that is allowed, and it is mirrored rather than
hot-linked because SourceForge's URLs are version-stamped and rot. The database is passed
with `-d`, but **only when this package fetched it**: a system ASTAP already knows where its
own database is, and overriding that would be presumptuous. It does its own star
detection, handles Alt-Az field rotation, and takes about a second per frame. This is the
right default: no rate limit, no network failure mode, no API key. Implementation notes:
- solves on a temporary 2-D copy of the **green plane**, with the header pointing and field
  of view as hints (a large robustness and speed win);
- runs a **two-pass cascade** — a fast downsampled solve (`-z 2`), then a full-resolution
  retry (`-z 0`) that recovers shallow stacks whose faint stars are lost at `-z 2`;
- solves in place with `-update`, then the clean WCS header is written to our standard
  sidecar, so it is a drop-in replacement for the astrometry.net path.

**`"nova"`** — astrometry.net's web service, fed **our own SEP source list** rather than the
image (much faster, and avoids re-detecting stars we already have), with the header pointing
and pixel scale as hints. Kept as a fallback for frames ASTAP cannot solve, and it is what
the historical MW Cam results were built with. It intermittently drops connections under
load, so only the network submit is retried, with backoff — never the extraction, which is
not a transient failure.

> **astrometry.net throttles above ~4 concurrent solves.** `pipeline.solve_all` therefore
> switches to a thread pool capped at 4 for `solver="nova"`, rather than one worker per core.
> More concurrency just produces dropped connections.

**`"lift"`** — for a frame that already carries a trustworthy solution. CrowdSky
plate-solves server-side and sets `PLTSOLVD = T`; `astrometry.lift` copies that header WCS
into a sidecar without solving. It **raises** if `PLTSOLVD` is absent, because lifting the
*Seestar's* on-board WCS would look like success and silently reintroduce the ~1 arcmin
errors.

## Degenerate solutions are never cached

Both solvers can return a non-celestial WCS. Those are rejected rather than written, and a
cached solution that turns out non-celestial is treated as absent, so a bad cache re-solves
transparently instead of poisoning the photometry.

`astrometry.has_wcs(path)` checks for a usable sidecar without loading the frame — the batch
runner uses it so that skipping a few thousand already-solved frames doesn't cost a few
thousand FITS reads.

## Pixel scale

`astrometry.pixel_scale(frame)` prefers a *measured* plate-solved scale
(`PIXSCALE`/`PLTSCALE`, set by CrowdSky) and falls back to the nominal
`206.265 × XPIXSZ / FOCALLEN`. The two agree to ~1% on the S50 (2.39 vs 2.37 arcsec/px).

## The reference catalogue

A Seestar dithers through a night, so a dataset covers more sky than any single FOV. One
**oversized catalogue is built once** and cached as an ECSV; every frame then subsets it to
its own footprint (`catalogs.sources_in_frame`). No per-frame catalogue queries — which is
what makes the per-frame work pure CPU and trivially parallel.

Defaults: a `2 × 1.5°` box, **one** cone, `phot_g_mean_mag < 17`.

- **Keep `n_tiles=1` unless the box is genuinely large.** Gaia TAP is unreliable under
  concurrent jobs; a single cone is both faster and safer.
- Tile centres are computed in a `SkyOffsetFrame`, not by adding degrees to RA — MW Cam sits
  at Dec +81, where naive RA steps collapse.
- The TAP endpoint intermittently truncates an async result (`IncompleteRead`), so the whole
  build is retried. **The cache is only written on full success**: a partial mosaic silently
  missing a tile is much worse than a failed build, because every frame in that region would
  lose its calibration stars with no error raised.

### Two backends, one cache

`Project.catalogue()` builds that ECSV either from the TAP query above or from an offline
copy of Gaia (`gaiadb`), and writes the same file either way, so nothing downstream can tell
which ran. `catalogue_backend="auto"` (the default) uses the offline copy when it covers the
field — installing the download is the only step needed to stop using the network.

The offline copy exists because **the TAP query is the least reliable thing the package
does**, and because a catalogue that is already on disk is what makes `solver="local"`
possible at all: a per-field query cannot solve the frame it is needed for.

Its shape, and the reasoning:

- **Rows are GSPC sources only** (those with synthetic photometry, G < 17.65), cut at
  V < 17.5. Restricting to sources that *have* synthetic V means every row is a usable
  calibrator with nothing masked in the columns that matter. Cutting on V rather than G
  costs only the ~10–15% of rows red enough for V to fall below G, and buys a faint limit
  that means something in the science band. About 195M rows.
- **Columns are what the package reads plus what it is about to want**, and no more: every
  extra float32 is ~0.5 GB over the full sky. Beyond the TAP set that means proper motions,
  `c_star` and `ipd_frac_multi_peak` (blend indicators — at 2.4 arcsec/px with a ~10 arcsec
  FWHM, blended comparisons are the norm), `ruwe` and `non_single_star` (unresolved
  binaries), `teff_gspphot`, and `v_jkc_flag`, which is Gaia's own statement that the
  synthetic V is inside its validated range.
- **Partitioned into 12288 HEALPix level-5 tiles**, so a partial copy is worth having: three
  observing fields cost tens of MB rather than the whole set. The partition key is free —
  a Gaia `source_id` already encodes its level-12 nested index, so `source_id >> 49` *is*
  the tile.

> **A cone query must read more sky than it asks for.** Completeness requires every tile
> whose *centre* lies within `radius + 1.91°`, because a tile is not a circle and its
> farthest corner is that far from its centre. Anything tighter silently drops sources
> instead of raising — and `astropy_healpix.cone_search_lonlat` is tighter, which is why it
> is not used. The exact cut afterwards removes the surplus, so the cost is I/O and never
> accuracy: measured over-read for a 1.5° cone is 12.4× the ideal at level 4, 5.1× at level
> 5, 2.7× at level 6. Level 5 is the compromise — level 6 would quadruple the file count to
> 49152 and shrink the parts below the size Parquet is efficient at.

### Proper motion

The TAP query does not fetch proper motions, so it cannot correct for them. Gaia DR3
positions are J2016.0; by 2026 a star with a 100 mas/yr proper motion has moved **1 arcsec**,
half the default 2 arcsec match tolerance, and the fastest movers are lost entirely. The
offline catalogue carries `pmra`/`pmdec`, and `Project(epoch=2026.4)` propagates positions
to the observing season. Left unset, positions stay at the Gaia epoch — today's behaviour.

### Why Gaia synthetic V

The query left-joins `gaiadr3.synthetic_photometry_gspc` for **synthetic
Johnson-Kron-Cousins V** (`v_jkc_mag`), plus `b_jkc_mag`/`r_jkc_mag` for a synthetic B−R
colour, and carries `phot_variable_flag` for vetting variables.

Gaia synthetic V beats APASS here: it is homogeneous across the sky, and its colour
information is what makes a per-frame colour term fittable at all. The green Bayer channel
through the Seestar's IRCUT filter is a good but inexact match to Johnson V; the residual
mismatch is exactly what the colour term absorbs.

Sources without synthetic photometry come through with those columns **masked**, so a
non-masked `v_jkc_mag` always means a usable calibrator.

### The zero-point fit

    V = m_inst + ZP + k·(B − R − colour0)

`ZP` is quoted at the field's median colour `colour0`, which decorrelates it from the colour
term — otherwise the two trade off and per-frame zero points are not comparable. The fit is
restricted to V ∈ [10, 14] (above saturation, below the detection floor), excludes
catalogue-flagged variables, and sigma-clips.

`calibration._pick_colour` cascades **B−R (JKC) → BP−RP → none**, and with no colour at all
degrades to a sigma-clipped zero-point mean with `colour_term = 0`. That exists so the fit
keeps working across catalogue schemas, not to squeeze out accuracy — the colour term is
small enough that a coarser colour barely moves the zero point.

## The API key

`ASTROMETRY_KEY` is read from the environment and never stored in the repo. It is only
needed for `solver="nova"`; the default ASTAP path is entirely offline.
