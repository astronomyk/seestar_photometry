"""Saved diagnostic figure sets -- one call per pipeline stage.

:mod:`plots` draws individual figures; this module composes them into named sets and
writes them to disk, so "did every stage do the right thing?" is a question you answer
by looking rather than by reasoning about numbers.

Three sets, matching the three pipeline stages:

:func:`frame_report`
    One frame, in depth: PSF, aperture, zero point, residual structure, background,
    detections, astrometry. Run it on a handful of frames -- this is where you find out
    *why* a dataset is behaving oddly.
:func:`frames_report`
    A whole dataset from its frame table: zero point and conditions over time, quality
    distributions, depth relations, plus a one-page contact sheet.
:func:`lightcurve_report`
    A finished light curve: finder chart, per-comparison curves, ensemble zero point,
    noise floor, periodogram, phase fold.

Filenames are deterministic, so re-running overwrites rather than accumulating and any
particular figure can be asked for by name. Every function returns the list of paths
written.
"""

from pathlib import Path

import numpy as np

from . import _style, plots

#: Figures are written at this size unless a function overrides it.
_DPI = 130


def _save(fig, outdir, name, written):
    """Write one figure, close it, and record the path."""
    import matplotlib.pyplot as plt

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{name}.png"
    fig.savefig(path, dpi=_DPI, bbox_inches="tight", facecolor=_style.SURFACE)
    plt.close(fig)
    written.append(path)
    return path


def _panel(draw, outdir, name, written, **kwargs):
    """Draw a single-axes figure via ``draw`` and save it, surviving a bad panel.

    A diagnostic set is most valuable exactly when something is wrong, which is also
    when an individual panel is most likely to fail (an empty band, a missing column, a
    degenerate fit). One dead panel must not cost you the other twenty, so the failure
    is drawn into the figure and the set continues.
    """
    import matplotlib.pyplot as plt

    try:
        result = draw(**kwargs)
        ax = result[0] if isinstance(result, tuple) else result
        fig = ax.figure if hasattr(ax, "figure") else ax
        return _save(fig, outdir, name, written)
    except Exception as exc:
        with plt.rc_context(_style.rc()):
            fig, ax = plt.subplots(figsize=(5.2, 3.6))
        ax.text(0.5, 0.5, f"{name}\ncould not be drawn:\n{exc}", ha="center",
                va="center", transform=ax.transAxes, fontsize=8,
                color=_style.STATUS["critical"], wrap=True)
        ax.set_axis_off()
        return _save(fig, outdir, name, written)


# --- per-frame ------------------------------------------------------------------------

def frame_report(frame, extraction, cal, outdir, prefix=None):
    """The full per-frame panel set.

    Parameters
    ----------
    frame : SeestarFrame
    extraction : Extraction
        Cross-matched (``match_gaia`` already applied).
    cal : Calibration
        From :func:`calibration.fit_zeropoint` on the same extraction, so it still
        carries its fit arrays.
    outdir : path-like
    prefix : str, optional
        Filename prefix. Defaults to the frame's stem, so several frames' reports
        coexist in one directory.
    """
    from . import photometry

    _style.use_style()
    outdir = Path(outdir)
    prefix = prefix or Path(frame.path).stem
    written = []

    _panel(plots.curve_of_growth, outdir, f"{prefix}_cog", written,
           cogs=extraction.cogs, apertures=extraction.aperture,
           enclosed=None)
    _panel(plots.fwhm_bands, outdir, f"{prefix}_fwhm", written, fwhm=extraction.fwhm)
    _panel(plots.reference_vs_instrumental, outdir, f"{prefix}_zp_relation", written,
           cal=cal)
    _panel(plots.zeropoint_vs_colour, outdir, f"{prefix}_zp_colour", written, cal=cal)
    for against in ("v", "snr", "radius"):
        _panel(plots.residual_vs, outdir, f"{prefix}_residual_{against}", written,
               cal=cal, against=against)
    _panel(plots.residual_map, outdir, f"{prefix}_residual_map", written,
           cal=cal, shape=frame.shape)
    _panel(plots.mag_snr, outdir, f"{prefix}_mag_snr", written,
           sources=extraction.sources)
    _panel(plots.match_separation, outdir, f"{prefix}_match_sep", written,
           sources=extraction.sources)
    _panel(plots.instrumental_colour, outdir, f"{prefix}_colour_fidelity", written,
           sources=extraction.sources)
    _panel(plots.detection_overlay, outdir, f"{prefix}_detections", written,
           frame=frame, extraction=extraction, half=250)

    # The background triptych is a multi-axes figure, so it is drawn directly.
    try:
        bg = photometry.fit_background(extraction, band="G")
        axes = plots.background_panels(bg)
        _save(axes[0].figure, outdir, f"{prefix}_background", written)
    except Exception:
        pass

    return written


def sample_frame_reports(project, n=3, keys=None):
    """Per-frame reports for the first ``n`` solved frames of a project.

    Re-measures those frames in the calling process rather than reusing the batch
    results, so diagnostics can be produced at any time -- including for a dataset whose
    frame table was already fully cached, which is the usual case when you come back to
    ask why something looks wrong.
    """
    from . import astrometry, calibration, frames, photometry

    catalogue = project.catalogue()
    keys = keys if keys is not None else project.frames()
    written = []
    made = 0
    for key in keys:
        if made >= n:
            break
        path = project.source.path(key)
        try:
            frame = frames.load_frame(path)
            wcs = astrometry.load_wcs(frame)
            if wcs is None:
                continue
            ext = photometry.extract_sources(
                frame, thresh=project.thresh, enclosed=project.enclosed_characterise
            )
            ext.match_gaia(catalogue, wcs=wcs, tol_arcsec=project.match_tol_arcsec)
            cal = calibration.fit_zeropoint(
                ext.sources, band="G", mag_range=project.fit_mag_range
            )
        except Exception as exc:
            print(f"[report] skipped {Path(str(path)).name}: {exc!r}", flush=True)
            continue
        written += frame_report(frame, ext, cal, project.diagnostics_dir)
        made += 1
    print(f"[report] wrote {len(written)} per-frame figures for {made} frames",
          flush=True)
    return written


# --- dataset --------------------------------------------------------------------------

def frames_report(frames_table, outdir, project=None):
    """Dataset-level panels from a frame table, plus a one-page contact sheet."""
    import matplotlib.pyplot as plt

    _style.use_style()
    outdir = Path(outdir)
    written = []
    if not len(frames_table):
        return written

    _panel(plots.zeropoint_vs_time, outdir, "frames_zp_vs_time", written,
           frames=frames_table)
    _panel(plots.metric_histogram, outdir, "frames_rms_hist", written,
           frames=frames_table, column="rms", threshold=0.06)
    _panel(plots.metric_histogram, outdir, "frames_chi2_hist", written,
           frames=frames_table, column="chi2_red", log=True)
    _panel(plots.depth_vs_exptime, outdir, "frames_depth_vs_exptime", written,
           frames=frames_table)
    for driver in ("sky_sb", "fwhm_G", "airmass"):
        if driver in frames_table.colnames:
            _panel(plots.depth_vs_driver, outdir, f"frames_depth_vs_{driver}", written,
                   frames=frames_table, driver=driver)
    for column in ("fwhm_G", "airmass", "sky_sb"):
        if column in frames_table.colnames:
            _panel(plots.condition_vs_time, outdir, f"frames_{column}_vs_time", written,
                   frames=frames_table, column=column)
    _panel(plots.calibration_coverage, outdir, "frames_coverage", written,
           frames=frames_table)

    # Contact sheet: the whole dataset at a glance, in one file.
    try:
        with plt.rc_context(_style.rc()):
            fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.4))
        plots.zeropoint_vs_time(frames_table, ax=axes[0, 0])
        plots.metric_histogram(frames_table, "rms", threshold=0.06, ax=axes[0, 1])
        plots.depth_vs_exptime(frames_table, ax=axes[0, 2])
        plots.condition_vs_time(frames_table, "fwhm_G", ax=axes[1, 0])
        plots.condition_vs_time(frames_table, "sky_sb", ax=axes[1, 1])
        plots.calibration_coverage(frames_table, ax=axes[1, 2])
        title = "Dataset summary"
        if project is not None:
            title += f" -- {project.target.name}"
        fig.suptitle(f"{title}   ({len(frames_table)} frames)", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        _save(fig, outdir, "frames_summary", written)
    except Exception as exc:
        print(f"[report] contact sheet failed: {exc!r}", flush=True)

    print(f"[report] wrote {len(written)} dataset figures to {outdir}", flush=True)
    return written


def measurements_report(stars, measurements, outdir, project=None, band="G"):
    """Quick panels straight off the measurement tables, before choosing an ensemble.

    Deliberately not the light-curve report: this runs with no comparison selection at
    all, so it can flag a problem (no target, no usable frames, a collapsed source
    count) before you spend any thought on which stars to compare against.
    """
    import matplotlib.pyplot as plt

    _style.use_style()
    outdir = Path(outdir)
    written = []
    target = stars[np.asarray(stars["is_target"])]
    if not len(target):
        return written
    target_id = int(target["source_id"][0])

    _panel(plots.raw_flux, outdir, "measurements_raw_flux", written,
           measurements=measurements, source_id=target_id, band=band)

    # Sources measured per frame: a collapse means frames that should be dropped.
    try:
        green = measurements[np.asarray(measurements["band"]) == band]
        frame_names, counts = np.unique(np.asarray(green["frame"]), return_counts=True)
        with plt.rc_context(_style.rc()):
            fig, ax = plt.subplots(figsize=(6.6, 3.0))
        ax.plot(np.arange(len(counts)), counts, ls="none", marker="o", markersize=3.5,
                color=_style.CATEGORICAL[0])
        ax.set_xlabel("frame index (alphabetical)")
        ax.set_ylabel("sources measured")
        ax.set_title(f"Catalogue sources per frame ({len(frame_names)} frames)")
        _save(fig, outdir, "measurements_sources_per_frame", written)
    except Exception as exc:
        print(f"[report] sources-per-frame failed: {exc!r}", flush=True)

    print(f"[report] wrote {len(written)} measurement figures to {outdir}", flush=True)
    return written


# --- light curve ----------------------------------------------------------------------

def lightcurve_report(lc, stars, measurements, comps, outdir, project=None,
                      band="G", frame=None, wcs=None, period=None, cutout=None,
                      contamination=None):
    """The full light-curve panel set.

    Parameters
    ----------
    lc : Table
        From :func:`lightcurves.differential_lightcurve`.
    stars, measurements : Table
        The two light-curve tables.
    comps : Table
        The comparison ensemble that produced ``lc``.
    outdir : path-like
    frame, wcs : optional
        A reference frame and its WCS, for the finder chart. Skipped if absent.
    period : float, optional
        Skip the period search and fold on this instead -- for a target with a known
        period, or to test a specific candidate.
    cutout, contamination : dict, optional
        From :mod:`contamination`, for a target sitting on extended emission.
    """
    from . import lightcurves

    _style.use_style()
    outdir = Path(outdir)
    written = []
    target_id = int(lc.meta.get("target_id", 0))
    name = project.target.name if project is not None else str(target_id)

    _panel(plots.lightcurve, outdir, "lc_differential", written,
           lc=lc, title=f"Differential light curve -- {name} ({band})")
    _panel(plots.ensemble_zeropoint, outdir, "lc_ensemble_zp", written, lc=lc)
    _panel(plots.raw_flux, outdir, "lc_raw_flux", written,
           measurements=measurements, source_id=target_id, band=band)

    if frame is not None and wcs is not None:
        _panel(plots.finder_chart, outdir, "lc_finder", written,
               frame=frame, wcs=wcs, stars=stars, comps=comps, target_id=target_id)

    # The per-comparison grid is the key check, so it is worth its own computation.
    curves = {}
    try:
        curves = lightcurves.comparison_curves(measurements, comps, band=band)
        if curves:
            fig, _axes = plots.comparison_grid(curves, comps=comps)
            _save(fig, outdir, "lc_comparison_grid", written)
            _panel(plots.scatter_vs_magnitude, outdir, "lc_noise_floor", written,
                   curves=curves, comps=comps)
    except Exception as exc:
        print(f"[report] comparison diagnostics failed: {exc!r}", flush=True)

    if len(lc) >= 5:
        try:
            pg = lightcurves.periodogram(lc)
            _panel(plots.periodogram, outdir, "lc_periodogram", written, pg=pg)
            fold_period = period if period is not None else pg["best_period"]
            _panel(plots.phase_fold, outdir, "lc_phase_fold", written,
                   lc=lc, period=fold_period)
        except Exception as exc:
            print(f"[report] period analysis failed: {exc!r}", flush=True)

    if cutout is not None:
        _panel(plots.host_cutout, outdir, "lc_host_cutout", written,
               cut=cutout, contamination=contamination)

    print(f"[report] wrote {len(written)} light-curve figures to {outdir}", flush=True)
    if curves:
        worst = sorted(curves.items(), key=lambda kv: -plots._scatter_of(kv[1]))[:3]
        print("[report] noisiest comparisons: "
              + ", ".join(f"{sid} ({plots._scatter_of(c) * 1000:.0f} mmag)"
                          for sid, c in worst), flush=True)
    return written
