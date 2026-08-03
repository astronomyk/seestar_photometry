# Astrometry and the reference catalogue

## The on-board WCS is unusable

Every Seestar frame arrives with a WCS in its header. It is **not** photometric-grade:
positions are wrong by tens of pixels (~1 arcmin), which drops the catalogue match rate to
a few percent and makes forced photometry meaningless — the aperture lands on empty sky.

So every frame is re-solved, and the solution is cached as a `.wcs` sidecar (a header-only
FITS) **next to the FITS file**. That location is deliberate and is the one exception to
"derived products go in `work_dir`": a solve is expensive, is a property of the frame rather
than of any analysis, and is reused by every project that touches that frame.

A good solve gives a **median cross-match separation of ~0.8 arcsec**. The
`*_match_sep` diagnostic panel is how you check: a distribution filling the match tolerance,
or pressed against it, means the solve is wrong even though it "succeeded".

## Two solvers, one sidecar

`astrometry.solve(frame, solver=...)`:

**`"astap"` (default)** — a local, fully offline blind plate solver. It does its own star
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
