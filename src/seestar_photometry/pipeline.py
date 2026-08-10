"""Batch stages: solve every frame, characterise every frame, measure every source.

Three stages, all driven from a :class:`project.Project`, all sharing one runner:

1. :func:`solve_all` -- a ``.wcs`` sidecar per frame.
2. :func:`build_frame_table` -- one quality row per frame, into ``frames.ecsv``.
3. :func:`build_measurements` -- forced photometry of every catalogue source in every
   frame, into ``stars.ecsv`` + ``measurements.ecsv``.

They must be run in that order: stages 2 and 3 read the cached WCS and never solve.
That separation is deliberate. Solving is network- or subprocess-bound and its cache
is per-frame and reusable; characterisation and measurement are pure CPU once the
WCS and catalogue are cached, so they parallelise cleanly across cores and can be
re-run freely as choices change.

Everything here is **resumable and idempotent**. Frames already in the output table
(or already carrying a sidecar) are skipped, results are checkpointed to disk
periodically, and a single bad frame is recorded with a status rather than sinking the
batch. Killing a run and restarting it is a normal operation, not a recovery
procedure -- which matters when a stage takes hours over thousands of frames.
"""

import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

#: Per-worker catalogue cache. Loading the ECSV once per process instead of once per
#: frame is the difference between seconds and minutes over a large dataset.
_CATALOGUE = None

#: Outcome of processing one frame.
#:
#: ``ok``          measured successfully
#: ``cached``      already done in a previous run, skipped
#: ``no_wcs``      no cached WCS sidecar -- run :func:`solve_all` first
#: ``load_error``  the FITS could not be read (truncated or malformed download)
#: ``failed``      raised during measurement; the message is logged
STATUSES = ("ok", "cached", "no_wcs", "load_error", "failed")


def _init_worker(catalogue_path):
    """Load the reference catalogue once per worker process, and quiet astroquery."""
    global _CATALOGUE
    from . import catalogs

    if catalogue_path is not None:
        _CATALOGUE = catalogs.load_catalogue(catalogue_path)
    try:
        from astroquery import log

        log.setLevel("ERROR")
    except Exception:
        pass


def _log(project, message):
    """Append a timestamped line to the project's error log, and echo it."""
    line = f"{time.strftime('%H:%M:%S')} {message}"
    print(line, flush=True)
    project.work_dir.mkdir(parents=True, exist_ok=True)
    with open(project.log_path, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def _run(project, work_fn, keys, workers=None, pool="process", catalogue_path=None,
         on_result=None, checkpoint=None, label="run", quiet=False):
    """Map ``work_fn`` over ``keys``, tolerating per-frame failure.

    ``work_fn(key, project)`` must return ``(key, status, payload)`` and must not
    raise -- failures are its own to catch and report, because a raise inside a worker
    loses the frame identity that makes the message actionable.

    ``pool`` is ``"process"`` for CPU-bound work or ``"thread"`` for network-bound
    work (astrometry.net wants a handful of concurrent submissions, not one per core).
    ``on_result(key, payload)`` is called in the parent for each ``ok``, and
    ``checkpoint(n_done)`` every ``checkpoint_every`` results.

    Returns ``(counts, problems)``.
    """
    workers = workers or (os.cpu_count() or 4)
    counts = dict.fromkeys(STATUSES, 0)
    problems = []
    # A single worker gains nothing from a separate process and pays the full spawn
    # cost -- which on Windows means re-importing astropy and scipy before any work
    # starts. Threads are strictly better there.
    if workers == 1:
        pool = "thread"
    executor = ProcessPoolExecutor if pool == "process" else ThreadPoolExecutor
    kwargs = {"max_workers": workers}
    if pool == "process":
        kwargs.update(initializer=_init_worker, initargs=(catalogue_path,))
    else:
        _init_worker(catalogue_path)  # threads share the parent's globals

    t0 = time.time()
    try:
        counts, problems = _drain(
            executor, kwargs, work_fn, keys, project, counts, problems,
            on_result, checkpoint, label, quiet,
        )
    except BrokenProcessPool as exc:
        # On Windows (spawn), a driver that calls a stage at module level instead of
        # inside `if __name__ == "__main__":` re-imports itself in every worker and the
        # pool dies with an opaque message. Say what to do about it.
        raise RuntimeError(
            "the worker pool died on startup. On Windows this almost always means the "
            "calling script runs the pipeline at import time -- wrap the call in\n"
            "    if __name__ == '__main__':\n        main()\n"
            "or pass workers=1 to run single-threaded."
        ) from exc

    if not quiet:
        summary = " ".join(f"{k} {v}" for k, v in counts.items() if v)
        print(f"[{label}] {summary} in {(time.time() - t0) / 60:.1f} min", flush=True)
    return counts, problems


def _drain(executor, kwargs, work_fn, keys, project, counts, problems,
           on_result, checkpoint, label, quiet):
    """Submit every key and consume results as they complete."""
    with executor(**kwargs) as ex:
        futures = {ex.submit(work_fn, key, project): key for key in keys}
        for done, future in enumerate(as_completed(futures), 1):
            key, status, payload = future.result()
            counts[status] = counts.get(status, 0) + 1
            if status == "ok" and on_result is not None:
                on_result(key, payload)
            elif status in ("failed", "load_error"):
                problems.append((status, key, payload))
                _log(project, f"{status}: {Path(str(key)).name} -- {payload}")
            if checkpoint is not None:
                checkpoint(done)
            if not quiet and done % 25 == 0:
                print(f"[{label} {done}/{len(futures)}] "
                      f"{ {k: v for k, v in counts.items() if v} }", flush=True)
    return counts, problems


# --- stage 1: WCS ---------------------------------------------------------------------

def _solve_one(key, project):
    """Solve and cache one frame's WCS. Returns ``(key, status, message)``."""
    from . import astrometry, frames

    path = project.source.path(key)
    try:
        if astrometry.has_wcs(path):
            return key, "cached", None
    except Exception:
        pass  # unreadable sidecar: fall through and re-solve

    try:
        frame = frames.load_frame(path)
    except Exception as exc:
        return key, "load_error", repr(exc)

    # Only the network/subprocess solve is retried; a failed extraction is not
    # transient, and re-running SEP on every attempt would multiply the cost of a
    # frame that was never going to solve. A local solve is deterministic, so retrying
    # it is pure waste -- it would fail identically.
    attempts = 4 if project.solver == "nova" else 1 if project.solver == "local" else 2
    last = None
    for attempt in range(attempts):
        try:
            astrometry.solve(
                frame, solver=project.solver, api_key=project.api_key,
                astap_exe=project.astap_exe, thresh=project.thresh,
                catalogue=_CATALOGUE,
            )
            return key, "solved", None
        except Exception as exc:
            last = repr(exc)
            if attempt < attempts - 1:
                time.sleep(5 * (attempt + 1))
    return key, "failed", last


def solve_all(project, workers=None, force=False, limit=None):
    """Solve and cache a WCS sidecar for every frame that lacks one.

    Idempotent: already-solved frames are skipped, so re-running fills gaps. The
    default ASTAP solver is a local subprocess and scales with cores; ``"nova"``
    switches to a small thread pool because astrometry.net accepts only a few
    concurrent jobs per account and answers the rest with dropped connections.

    Returns ``(counts, problems)``.
    """
    from . import astrometry

    project.ensure_dirs()
    keys = project.frames()
    if not force:
        keys = [k for k in keys if not astrometry.has_wcs(project.source.path(k))]
    if limit is not None:
        keys = keys[:limit]
    print(f"[solve] {len(keys)} frames need solving (solver={project.solver})", flush=True)
    if not keys:
        return dict.fromkeys(STATUSES, 0), []

    # The local solver pairs against the reference catalogue, so unlike the other
    # backends this stage needs it -- built once here in the parent, then loaded once
    # per worker like the later stages do.
    catalogue_path = None
    if project.solver == "local":
        catalogue_path = project.catalogue_path
        if not catalogue_path.exists():
            project.catalogue()

    if project.solver == "nova":
        pool, workers = "thread", min(workers or 4, 4)
    else:
        pool, workers = "process", workers or max((os.cpu_count() or 4) - 2, 1)
    counts, problems = _run(
        project, _solve_one, keys, workers=workers, pool=pool, label="solve",
        catalogue_path=catalogue_path,
    )
    return counts, problems


# --- stage 2: per-frame characterisation ----------------------------------------------

def _frame_mask(project, frame, wcs):
    """The project's per-frame exclusion mask, or ``None``.

    ``Project.mask`` is a callable so the region can be specified on the sky once and
    resolved to pixels per frame -- dithering and Alt-Az rotation move it. A mask that
    raises is treated as fatal rather than skipped: quietly measuring an unmasked frame
    alongside masked ones would mix two aperture-sizing regimes in one table.
    """
    if getattr(project, "mask", None) is None:
        return None
    return project.mask(frame, wcs)


def _frame_row(key, project):
    """Characterise one frame into a quality row. Returns ``(key, status, row)``."""
    from . import astrometry, calibration, frames, photometry, quality

    path = project.source.path(key)
    try:
        frame = frames.load_frame(path)
    except Exception as exc:
        return key, "load_error", repr(exc)
    wcs = astrometry.load_wcs(frame)
    if wcs is None:
        return key, "no_wcs", None
    try:
        ext = photometry.extract_sources(
            frame, thresh=project.thresh, enclosed=project.enclosed_characterise,
            mask=_frame_mask(project, frame, wcs), isolation=project.isolation,
        )
        ext.match_gaia(_CATALOGUE, wcs=wcs, tol_arcsec=project.match_tol_arcsec)
        cal = calibration.fit_zeropoint(
            ext.sources, band="G", mag_range=project.fit_mag_range
        )
        row = quality.frame_quality(ext, cal, provenance=project.provenance)
        row.update(astrometry.match_quality(ext.band("G"), project.match_tol_arcsec))
        row.update(quality.onboard_quality(frame))
        row["solver"] = project.solver
        row["enclosed"] = project.enclosed_characterise
        return key, "ok", row
    except Exception as exc:
        return key, "failed", f"{exc!r}\n{traceback.format_exc()}"


def build_frame_table(project, workers=None, force=False, limit=None,
                      checkpoint_every=25, diagnostics=False):
    """Characterise every solved frame into ``frames.ecsv``.

    Reads the cached WCS and never solves -- run :func:`solve_all` first. Rows already
    in the table are skipped unless ``force``, and the table is rewritten every
    ``checkpoint_every`` frames so an interrupted run loses at most that many.

    ``diagnostics`` saves figures into ``work_dir/diagnostics/``: ``True`` for the
    dataset-level panels only, an ``int`` to also save the per-frame panel set for
    that many frames. Per-frame panels re-measure their frames in the parent process,
    so they cost a little extra time but work even when every row was already cached.

    Returns the full frame table.
    """
    project.ensure_dirs()
    catalogue_path = project.catalogue_path
    if not catalogue_path.exists():
        project.catalogue()  # build it once, in the parent, before forking

    existing = None
    done = set()
    if project.frames_path.exists() and not force:
        existing = Table.read(project.frames_path)
        done = {os.path.normcase(str(p)) for p in existing["path"]}

    keys = [k for k in project.frames()
            if os.path.normcase(str(project.source.path(k))) not in done]
    if limit is not None:
        keys = keys[:limit]
    print(f"[frames] {len(keys)} to measure ({len(done)} already done)", flush=True)

    rows = []
    def _write(_n_done=None):
        if not rows:
            return
        new = Table(rows)
        table = vstack([existing, new], metadata_conflicts="silent") \
            if existing is not None else new
        table.write(project.frames_path, format="ascii.ecsv", overwrite=True)

    if keys:
        _run(
            project, _frame_row, keys, workers=workers, pool="process",
            catalogue_path=catalogue_path, label="frames",
            on_result=lambda k, row: rows.append(row),
            checkpoint=lambda n: _write() if n % checkpoint_every == 0 else None,
        )
        _write()

    table = Table.read(project.frames_path) if project.frames_path.exists() else Table()
    if len(table):
        table.sort(["unit", "date_obs"])
        table.write(project.frames_path, format="ascii.ecsv", overwrite=True)
        _report_frame_summary(table)

    if diagnostics:
        from . import report

        report.frames_report(table, project.diagnostics_dir, project=project)
        n_frames = diagnostics if isinstance(diagnostics, int) and diagnostics is not True else 0
        if n_frames:
            report.sample_frame_reports(project, n=n_frames)
    return table


def _report_frame_summary(table):
    """Print the numbers you actually want after a characterisation run."""
    rms = np.asarray(table["rms"], dtype=float)
    p10, p50, p90 = np.nanpercentile(rms, [10, 50, 90])
    print(f"[frames] {len(table)} rows | calibration rms 10/50/90 pct = "
          f"{p10:.3f} / {p50:.3f} / {p90:.3f} | "
          f"rms < 0.06: {int((rms < 0.06).sum())}/{len(table)}", flush=True)


# --- stage 3: forced photometry of every source ---------------------------------------

def _measure_sources(key, project):
    """Forced photometry + calibration for one frame. Returns a measurement Table."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    from . import astrometry, calibration, catalogs, frames, lightcurves, photometry

    path = project.source.path(key)
    try:
        frame = frames.load_frame(path)
    except Exception as exc:
        return key, "load_error", repr(exc)
    wcs = astrometry.load_wcs(frame)
    if wcs is None:
        return key, "no_wcs", None
    try:
        in_frame = catalogs.sources_in_frame(_CATALOGUE, wcs, frame)
        forced = photometry.forced_photometry(
            frame, in_frame["ra"], in_frame["dec"], wcs,
            source_id=in_frame["source_id"], thresh=project.thresh,
            enclosed=project.enclosed_lightcurve,
            mask=_frame_mask(project, frame, wcs), isolation=project.isolation,
        )
        # The forced positions *are* catalogue sources, so their catalogue photometry
        # is carried straight through rather than re-matched. `forced` is three
        # per-band blocks, each in `in_frame` row order.
        n_band = len(photometry.BANDS)
        for col in ("v_jkc_mag", "b_jkc_mag", "r_jkc_mag", "phot_variable_flag"):
            if col in in_frame.colnames:
                forced[col] = np.ma.concatenate([in_frame[col]] * n_band)

        cal = calibration.fit_zeropoint(
            forced, band="G", mag_range=project.fit_mag_range
        )

        target = SkyCoord(project.target.ra * u.deg, project.target.dec * u.deg)
        meta = frames.frame_metadata(frame)
        times = lightcurves.frame_times(frame.header, target, meta["total_exptime"])

        m_inst = photometry.instrumental_mag(forced["flux"])
        mag_err = photometry.mag_error(forced["snr"])
        # Only the green band is calibrated onto reference V; R and B keep their
        # instrumental magnitudes, which is all the colour diagnostics need.
        colour = (np.asarray(forced["b_jkc_mag"], dtype=float)
                  - np.asarray(forced["r_jkc_mag"], dtype=float)) \
            if "b_jkc_mag" in forced.colnames else None
        is_green = np.asarray(forced["band"]) == "G"
        mag = np.where(
            is_green, calibration.apply_calibration(m_inst, cal, colour), np.nan
        )

        n = len(forced)
        t = Table()
        t["source_id"] = forced["source_id"]
        t["frame"] = np.full(n, Path(str(path)).name)
        t["band"] = forced["band"]
        t["mjd_obs"] = np.full(n, times["mjd_obs"])
        t["mjd_mid"] = np.full(n, times["mjd_mid"])
        t["bjd_tdb"] = np.full(n, times["bjd_tdb"])
        t["airmass"] = np.full(n, meta["airmass"])
        t["x"] = forced["x"]
        t["y"] = forced["y"]
        t["flux"] = forced["flux"]
        t["fluxerr"] = forced["fluxerr"]
        t["snr"] = forced["snr"]
        t["m_inst"] = m_inst
        t["mag"] = mag
        t["mag_err"] = mag_err
        t["flag"] = forced["flag"]
        t["on_chip"] = forced["on_chip"]
        t["max_pix_value"] = forced["max_pix_value"]
        t["zeropoint"] = np.full(n, cal.zeropoint)
        t["zp_rms"] = np.full(n, cal.rms)
        return key, "ok", t
    except Exception as exc:
        return key, "failed", f"{exc!r}\n{traceback.format_exc()}"


def build_measurements(project, workers=None, limit=None, diagnostics=False):
    """Forced-photometer every catalogue source in every solved frame.

    Writes ``measurements.ecsv`` (one row per source, frame and band) and
    ``stars.ecsv`` (one row per source). Unlike :func:`build_frame_table` this is not
    incremental -- the long table is rebuilt in one pass, because it is cheap relative
    to solving and because a partially-rebuilt light curve mixing two aperture
    choices would be a subtle, invisible error.

    Returns ``(stars, measurements)``.
    """
    project.ensure_dirs()
    catalogue_path = project.catalogue_path
    if not catalogue_path.exists():
        project.catalogue()

    keys = project.frames()
    if limit is not None:
        keys = keys[:limit]
    print(f"[measure] {len(keys)} frames", flush=True)

    parts = []
    _run(
        project, _measure_sources, keys, workers=workers, pool="process",
        catalogue_path=catalogue_path, label="measure",
        on_result=lambda k, t: parts.append(t),
    )
    if not parts:
        raise RuntimeError(
            "no frames were measured; run solve_all() first and check errors.log"
        )

    from . import catalogs, lightcurves

    measurements = vstack(parts, metadata_conflicts="silent")
    measurements.sort(["source_id", "band", "mjd_mid"])
    catalogue = catalogs.load_catalogue(catalogue_path)
    stars = lightcurves.build_stars(
        catalogue, np.unique(measurements["source_id"]), project.target.radec,
        target_id=project.target.source_id,
    )
    measurements.write(project.measurements_path, format="ascii.ecsv", overwrite=True)
    stars.write(project.stars_path, format="ascii.ecsv", overwrite=True)

    print(f"[measure] {len(measurements)} measurements of {len(stars)} sources "
          f"-> {project.measurements_path}", flush=True)
    target = stars[np.asarray(stars["is_target"])]
    if len(target):
        print(f"[measure] target source_id={int(target['source_id'][0])} "
              f"V={float(target['v_jkc_mag'][0]):.2f} "
              f"sep={float(target['sep_target_arcmin'][0]):.2f}'", flush=True)

    if diagnostics:
        from . import report

        report.measurements_report(
            stars, measurements, project.diagnostics_dir, project=project
        )
    return stars, measurements


def load_tables(project):
    """Read back whichever of the three output tables exist.

    Returns ``(frames, stars, measurements)`` with ``None`` for any not yet built.
    """
    read = lambda p: Table.read(p) if Path(p).exists() else None
    return (read(project.frames_path), read(project.stars_path),
            read(project.measurements_path))
