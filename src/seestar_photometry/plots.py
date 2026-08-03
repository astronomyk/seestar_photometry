"""Individual diagnostic figures. One function per figure, no file I/O.

Each function takes the data it needs plus an optional ``ax``, draws into it and
returns the axes -- so the same figure composes into a notebook, into a multi-panel
sheet, or into a saved PNG (see :mod:`report`) without changing.

Matplotlib is imported lazily inside each function, so ``import seestar_photometry``
works in a core install with no plotting dependency.

Styling, and in particular the per-band colour/marker/linestyle convention and why
the secondary encoding is mandatory, lives in :mod:`_style`.
"""

import numpy as np

from . import _style
from ._style import new_axes

#: Alpha for dense point clouds. A field is hundreds to thousands of stars, which is a
#: density distribution rather than a set of labelled series -- so these are drawn as
#: small translucent points instead of the large opaque marks a categorical series
#: gets. Overplotting, not mark size, is the thing to solve here.
_CLOUD = {"s": 9, "alpha": 0.45, "linewidths": 0}


# --- per-frame: PSF and aperture ------------------------------------------------------

def curve_of_growth(cogs, apertures=None, enclosed=None, ax=None):
    """Curve of growth per band, with the chosen aperture radius marked.

    Confirms the aperture sizing is working: each curve should rise smoothly and
    flatten, and the marked radius should sit on the shoulder. A curve that keeps
    climbing to the last radius means the "total" normalisation is contaminated by a
    neighbour, and the aperture will have been sized too large for the whole frame.
    """
    from ._style import BAND_COLOR

    fig, ax = new_axes(ax)
    for i, band in enumerate(("R", "G", "B")):
        cog = cogs[i] if cogs is not None and i < len(cogs) else None
        if cog is None:
            continue
        radius = np.asarray(cog["radius"], dtype=float)
        frac = np.asarray(cog["flux_frac"], dtype=float)
        std = np.asarray(cog["flux_frac_std"], dtype=float)
        ax.plot(radius, frac, ls=_style.BAND_STYLE[band], **_style.band_kw(band),
                markersize=3.5)
        ax.fill_between(radius, frac - std, frac + std,
                        color=BAND_COLOR[band], alpha=0.10, lw=0)
        if apertures is not None and np.isfinite(apertures[i]):
            ax.axvline(apertures[i], color=BAND_COLOR[band], lw=0.9, ls="-",
                       alpha=0.55, zorder=1)
    if enclosed is not None:
        _style.reference_line(ax, enclosed, axis="y", label=f"{enclosed:.0%} enclosed")
    ax.set_xlabel("aperture radius (pixels)")
    ax.set_ylabel("enclosed flux fraction")
    ax.set_title("Curve of growth, per band")
    ax.set_ylim(0, 1.05)
    _style.legend(ax, title="band")
    if apertures is not None:
        _style.annotate(
            ax,
            "aperture (px)\n" + "\n".join(
                f"{b}  {apertures[i]:.2f}" for i, b in enumerate(("R", "G", "B"))
                if np.isfinite(apertures[i])
            ),
            loc="lower right",
        )
    return ax


def fwhm_bands(fwhm, ax=None):
    """Per-band FWHM and the R/G and B/G ratios.

    Tests the chromatic-PSF premise the per-band aperture rests on: the Seestar focuses
    best near green, so R and B should come out broader (roughly R/G 1.04-1.20,
    B/G 1.06-1.30). Ratios near 1.00 in every frame would mean the per-band sizing is
    buying nothing; ratios far above the range point at bad focus or trailing.
    """
    fig, ax = new_axes(ax, figsize=(4.2, 3.4))
    bands = ("R", "G", "B")
    values = [float(fwhm[i]) for i in range(3)]
    colors = [_style.BAND_COLOR[b] for b in bands]
    # 2px surface gap between adjacent bars keeps the fills from touching.
    ax.bar(bands, values, color=colors, width=0.62, edgecolor=_style.SURFACE, lw=2)
    for i, v in enumerate(values):
        ax.text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=8,
                color=_style.INK_SECONDARY)
    ax.set_ylabel("FWHM (pixels)")
    ax.set_title("PSF width per band")
    ax.grid(axis="x", visible=False)
    g = values[1]
    if g > 0:
        _style.annotate(ax, f"R/G  {values[0] / g:.3f}\nB/G  {values[2] / g:.3f}")
    return ax


# --- per-frame: zero point ------------------------------------------------------------

def _fit_arrays(cal):
    """The fit arrays a calibration carries, or raise a useful error."""
    if not cal.fit:
        raise ValueError(
            "this Calibration carries no fit arrays; it was probably read back from "
            "a table rather than produced by fit_zeropoint()"
        )
    return cal.fit


def reference_vs_instrumental(cal, ax=None):
    """Reference V against instrumental magnitude -- the zero-point relation.

    Should be a straight line of slope 1. The fit window is shaded; points outside it
    show whether the window is in the right place. A flattening at the bright end is
    saturation, and a fan-out at the faint end is the noise floor -- if either intrudes
    into the shaded window, move the window rather than trusting the fit.
    """
    fit = _fit_arrays(cal)
    fig, ax = new_axes(ax)
    keep = fit["keep"]
    ax.scatter(fit["m_inst"][~keep], fit["v"][~keep], color=_style.INK_MUTED,
               label="sigma-clipped", **_CLOUD)
    ax.scatter(fit["m_inst"][keep], fit["v"][keep], color=_style.BAND_COLOR[cal.band],
               label="used in fit", **_CLOUD)
    if cal.mag_range:
        ax.axhspan(cal.mag_range[0], cal.mag_range[1], color=_style.MODEL,
                   alpha=0.07, lw=0, zorder=0)
        ax.text(0.02, cal.mag_range[0], " fit window",
                transform=ax.get_yaxis_transform(), fontsize=7,
                color=_style.REFERENCE, va="bottom")
    x = np.linspace(np.nanmin(fit["m_inst"]), np.nanmax(fit["m_inst"]), 10)
    ax.plot(x, x + cal.zeropoint, color=_style.MODEL, lw=1.4,
            label=f"ZP = {cal.zeropoint:.3f}")
    ax.set_xlabel("instrumental magnitude  ($-2.5\\log_{10}$ flux)")
    ax.set_ylabel("reference V (mag)")
    ax.set_title(f"Zero-point relation, {cal.band} band")
    _style.mag_axis(ax, "y")
    _style.legend(ax, loc="lower left")
    _style.annotate(ax, f"n = {cal.n_stars}\nrms = {cal.rms:.3f} mag", loc="upper right")
    return ax


def zeropoint_vs_colour(cal, ax=None):
    """Zero-point residual against colour -- the colour-term fit.

    The residual ``V - m_inst`` is ``ZP + k*(colour - colour0)``, so this panel *is*
    the calibration. The fitted line's slope is the colour term, which exists because
    the green Bayer channel through IRCUT is not exactly Johnson V. A slope that
    changes a lot frame to frame means something is varying that a constant colour
    term can't absorb -- usually focus, via the chromatic PSF.
    """
    fit = _fit_arrays(cal)
    if fit["colour"] is None:
        fig, ax = new_axes(ax)
        ax.text(0.5, 0.5, "no colour available\n(zero point fit as a clipped mean)",
                ha="center", va="center", transform=ax.transAxes,
                color=_style.INK_SECONDARY, fontsize=9)
        ax.set_axis_off()
        return ax

    fig, ax = new_axes(ax)
    colour, resid, keep = fit["colour"], fit["residual"], fit["keep"]
    ax.scatter(colour[~keep], resid[~keep], color=_style.INK_MUTED,
               label="sigma-clipped", **_CLOUD)
    ax.scatter(colour[keep], resid[keep], color=_style.BAND_COLOR[cal.band],
               label="used in fit", **_CLOUD)
    x = np.linspace(np.nanmin(colour), np.nanmax(colour), 10)
    ax.plot(x, cal.zeropoint + cal.colour_term * (x - cal.colour0),
            color=_style.MODEL, lw=1.4, label=f"k = {cal.colour_term:+.3f}")
    _style.reference_line(ax, cal.colour0, axis="x", label="colour0")
    ax.set_xlabel(f"colour  {cal.colour_label}")
    ax.set_ylabel("$V - m_{inst}$  (mag)")
    ax.set_title("Colour term")
    _style.legend(ax, loc="lower left")
    _style.annotate(
        ax,
        f"ZP = {cal.zeropoint:.3f} ± {cal.zeropoint_err:.3f}\n"
        f"k  = {cal.colour_term:+.3f} ± {cal.colour_term_err:.3f}\n"
        f"rms = {cal.rms:.3f} mag",
        loc="upper right",
    )
    return ax


def residual_vs(cal, against="v", ax=None):
    """Zero-point residual against magnitude, SNR or radius from the frame centre.

    All three should be flat. Structure means a systematic the zero point is absorbing
    into its scatter:

    * against magnitude -- non-linearity, or saturation creeping into the fit window;
    * against SNR -- the noise model is wrong;
    * against radius -- vignetting or flat-field residual, which on an Alt-Az mount
      also shows up as field-rotation-dependent aperture loss.
    """
    fit = _fit_arrays(cal)
    keep = fit["keep"]
    resid = fit["residual"] - fit["model"]
    if against == "v":
        x, label = fit["v"], "reference V (mag)"
    elif against == "snr":
        x, label = fit["snr"], "SNR"
    elif against == "radius":
        cx, cy = np.nanmedian(fit["x"]), np.nanmedian(fit["y"])
        x = np.hypot(fit["x"] - cx, fit["y"] - cy)
        label = "radius from frame centre (pixels)"
    else:
        raise ValueError(f"unknown x-axis {against!r} (expected v, snr or radius)")

    fig, ax = new_axes(ax)
    ax.scatter(x[~keep], resid[~keep], color=_style.INK_MUTED, **_CLOUD)
    ax.scatter(x[keep], resid[keep], color=_style.BAND_COLOR[cal.band], **_CLOUD)
    _style.reference_line(ax, 0.0, axis="y")
    if against == "snr":
        ax.set_xscale("log")
    ax.set_xlabel(label)
    ax.set_ylabel("residual (mag)")
    ax.set_title(f"Residual vs {against}")
    # A running median makes a weak trend visible that a point cloud hides.
    ok = np.isfinite(x) & np.isfinite(resid) & keep
    if ok.sum() > 20:
        order = np.argsort(x[ok])
        xs, ys = x[ok][order], resid[ok][order]
        n_bin = max(int(len(xs) / 20), 4)
        edges = np.array_split(np.arange(len(xs)), n_bin)
        bx = [np.median(xs[e]) for e in edges if len(e)]
        by = [np.median(ys[e]) for e in edges if len(e)]
        ax.plot(bx, by, color=_style.MODEL, lw=1.4, label="running median")
        ax.legend(loc="upper right")
    return ax


def residual_map(cal, shape=None, ax=None):
    """Zero-point residual laid out over the frame, as a diverging map.

    A flat field of mixed signs is what you want. Any coherent structure -- a gradient,
    a corner, a ring -- is a spatial systematic, and it will not cancel in differential
    photometry unless the comparisons happen to sit where the target does. This is the
    panel that justifies a proximity cut when choosing comparisons.
    """
    fit = _fit_arrays(cal)
    keep = fit["keep"]
    resid = (fit["residual"] - fit["model"])[keep]
    fig, ax = new_axes(ax, figsize=(4.6, 4.0))
    # Diverging scale: two hue poles about a neutral zero, symmetric limits so the
    # midpoint reads as "no residual".
    span = float(np.nanpercentile(np.abs(resid), 95)) or 0.05
    art = ax.scatter(fit["x"][keep], fit["y"][keep], c=resid, cmap="RdBu_r",
                     vmin=-span, vmax=span, s=16, linewidths=0)
    bar = ax.figure.colorbar(art, ax=ax, fraction=0.046, pad=0.03)
    bar.set_label("residual (mag)", fontsize=8, color=_style.INK_SECONDARY)
    bar.outline.set_visible(False)
    if shape is not None:
        ax.set_xlim(0, shape[1])
        ax.set_ylim(0, shape[0])
    ax.set_aspect("equal")
    ax.set_xlabel("x (pixels)")
    ax.set_ylabel("y (pixels)")
    ax.set_title("Residual across the frame")
    ax.grid(visible=False)
    return ax


# --- per-frame: detection and astrometry ----------------------------------------------

def mag_snr(sources, ax=None):
    """Instrumental magnitude against SNR, per band.

    Shows where the SNR > 5 sample boundary falls and how the three bands compare in
    depth. A second, offset streak of points is a sign of two source populations being
    measured differently -- the known signature of the Kron minimum-aperture fallback,
    for instance.
    """
    fig, ax = new_axes(ax)
    band_col = np.asarray(sources["band"])
    mag = -2.5 * np.log10(np.asarray(sources["flux"], dtype=float))
    snr = np.asarray(sources["snr"], dtype=float)
    for band in ("R", "G", "B"):
        m = band_col == band
        ax.scatter(mag[m], snr[m], color=_style.BAND_COLOR[band],
                   marker=_style.BAND_MARKER[band], label=band, **_CLOUD)
    _style.reference_line(ax, 5.0, axis="y", label="SNR 5")
    ax.set_yscale("log")
    ax.set_xlabel("instrumental magnitude  ($-2.5\\log_{10}$ flux)")
    ax.set_ylabel("SNR")
    ax.set_title("Depth per band")
    _style.legend(ax, title="band")
    return ax


def background_panels(bg, axes=None):
    """SEP background mesh, the fitted 2nd-order polynomial, and the residual.

    The polynomial should absorb the smooth gradient, leaving a residual that is
    structureless apart from halos around bright stars -- those genuinely belong to the
    sources, not the sky. Large-scale structure left in the residual means the sky is
    not describable by a quadratic, and the pedestal (hence the sky-brightness number)
    is being biased.
    """
    import matplotlib.pyplot as plt

    if axes is None:
        with plt.rc_context(_style.rc()):
            fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    lo, hi = np.nanpercentile(bg["background"], [1, 99])
    for ax, (image, title) in zip(
        axes,
        ((bg["background"], "SEP background mesh"),
         (bg["model"], "2nd-order polynomial"),
         (bg["residual"], "residual")),
    ):
        if title == "residual":
            span = float(np.nanpercentile(np.abs(image), 99)) or 1.0
            art = ax.imshow(image, origin="lower", cmap="RdBu_r",
                            vmin=-span, vmax=span)
        else:
            art = ax.imshow(image, origin="lower", cmap="Blues", vmin=lo, vmax=hi)
        bar = ax.figure.colorbar(art, ax=ax, fraction=0.046, pad=0.03)
        bar.outline.set_visible(False)
        bar.ax.tick_params(labelsize=7)
        ax.set_title(title)
        ax.grid(visible=False)
        ax.set_xticks([])
        ax.set_yticks([])
    _style.annotate(
        axes[0],
        f"pedestal {bg['pedestal']:.1f} ADU\nresid std {bg['resid_std']:.2f}",
    )
    return axes


def detection_overlay(frame, extraction, ax=None, matched_only=False, half=None):
    """The green plane with measured apertures drawn on, and catalogue matches ringed.

    The direct check that astrometry, detection and cross-match agree. Apertures should
    sit centred on stars, and matched sources should carry a ring; a systematic offset
    between rings and apertures is a bad WCS, which is far easier to see here than in
    any summary statistic.
    """
    fig, ax = new_axes(ax, figsize=(5.0, 5.0))
    green = np.asarray(frame.g, dtype=float)
    ny, nx = green.shape
    if half:
        cy, cx = ny // 2, nx // 2
        y0, y1 = max(cy - half, 0), min(cy + half, ny)
        x0, x1 = max(cx - half, 0), min(cx + half, nx)
    else:
        y0, y1, x0, x1 = 0, ny, 0, nx
    stamp = green[y0:y1, x0:x1]
    lo, hi = np.nanpercentile(stamp, [30, 99.5])
    ax.imshow(stamp, origin="lower", cmap="gray_r", vmin=lo, vmax=hi,
              extent=(x0, x1, y0, y1))

    g = extraction.band("G")
    radius = float(extraction.aperture[1]) if extraction.aperture is not None else 4.0
    matched = (~np.ma.getmaskarray(g["v_jkc_mag"])) if "v_jkc_mag" in g.colnames \
        else np.zeros(len(g), dtype=bool)
    x, y = np.asarray(g["x"], dtype=float), np.asarray(g["y"], dtype=float)
    inside = (x >= x0) & (x < x1) & (y >= y0) & (y < y1)
    show = inside & matched if matched_only else inside

    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Circle

    ax.add_collection(PatchCollection(
        [Circle((xi, yi), radius) for xi, yi in zip(x[show], y[show])],
        facecolors="none", edgecolors=_style.BAND_COLOR["G"], linewidths=0.8,
    ))
    ring = show & matched
    ax.add_collection(PatchCollection(
        [Circle((xi, yi), radius * 2.1) for xi, yi in zip(x[ring], y[ring])],
        facecolors="none", edgecolors=_style.CATEGORICAL[1], linewidths=0.7,
        linestyles="dotted",
    ))
    ax.set_title("Apertures (green) and catalogue matches (dotted)")
    ax.set_xlabel("x (pixels)")
    ax.set_ylabel("y (pixels)")
    ax.grid(visible=False)
    _style.annotate(
        ax, f"aperture {radius:.2f} px\n{int(ring.sum())} matched / {int(show.sum())} shown"
    )
    return ax


def match_separation(sources, tol_arcsec=2.0, ax=None):
    """Histogram of cross-match separations -- the practical test of the WCS solve.

    A good solve piles up below ~1 arcsec. A broad distribution filling the tolerance,
    or a peak pressed against it, means the solve is wrong even though it "succeeded" --
    the failure mode the unusable on-board WCS produces, and the reason every frame is
    re-solved.
    """
    fig, ax = new_axes(ax, figsize=(4.6, 3.4))
    sep = np.asarray(sources["sep_arcsec"], dtype=float)
    sep = sep[np.isfinite(sep) & (sep < tol_arcsec * 3)]
    ax.hist(sep, bins=40, color=_style.CATEGORICAL[0], edgecolor=_style.SURFACE, lw=0.6)
    _style.reference_line(ax, tol_arcsec, axis="x", label="match tolerance")
    ax.set_xlabel("separation to nearest catalogue source (arcsec)")
    ax.set_ylabel("sources")
    ax.set_title("Cross-match quality")
    if len(sep):
        _style.annotate(
            ax,
            f"median {np.median(sep):.2f}\"\n"
            f"90th pct {np.percentile(sep, 90):.2f}\"\n"
            f"n = {len(sep)}",
            loc="upper right",
        )
    return ax


def instrumental_colour(sources, ax=None):
    """Instrumental B-R against catalogue B-R, per matched star.

    Tests the per-band aperture directly: since every band encloses the same flux
    fraction, the aperture correction cancels in a colour, so instrumental colour
    should track catalogue colour with unit slope and constant offset. Curvature or a
    slope away from 1 means the bands are *not* enclosing matched fractions, which is
    exactly the bias a single shared radius would introduce.
    """
    fig, ax = new_axes(ax)
    band_col = np.asarray(sources["band"])
    need = {"b_jkc_mag", "r_jkc_mag"}
    if not need <= set(sources.colnames):
        ax.text(0.5, 0.5, "no catalogue colour available", ha="center", va="center",
                transform=ax.transAxes, color=_style.INK_SECONDARY)
        ax.set_axis_off()
        return ax

    # Pair R and B rows by position: the per-band blocks are in the same source order.
    r = sources[band_col == "R"]
    b = sources[band_col == "B"]
    n = min(len(r), len(b))
    inst = (-2.5 * np.log10(np.asarray(b["flux"][:n], dtype=float))
            + 2.5 * np.log10(np.asarray(r["flux"][:n], dtype=float)))
    cat = (np.asarray(b["b_jkc_mag"][:n], dtype=float)
           - np.asarray(b["r_jkc_mag"][:n], dtype=float))
    ok = np.isfinite(inst) & np.isfinite(cat)
    ax.scatter(cat[ok], inst[ok], color=_style.CATEGORICAL[0], **_CLOUD)
    if ok.sum() > 5:
        slope, intercept = np.polyfit(cat[ok], inst[ok], 1)
        x = np.linspace(np.nanmin(cat[ok]), np.nanmax(cat[ok]), 10)
        ax.plot(x, slope * x + intercept, color=_style.MODEL, lw=1.4,
                label=f"slope {slope:.3f}")
        ax.legend(loc="lower right")
    ax.set_xlabel("catalogue B-R (mag)")
    ax.set_ylabel("instrumental B-R (mag)")
    ax.set_title("Colour fidelity of the per-band aperture")
    return ax


# --- dataset level --------------------------------------------------------------------

def _time_axis(frames):
    """Relative x values and an axis label, from whatever time column exists."""
    from astropy.time import Time

    raw = None
    for col in ("mjd_mid", "mjd_obs"):
        if col in frames.colnames:
            raw = np.asarray(frames[col], dtype=float)
            break
    if raw is None:
        iso = [str(v).replace(" ", "T") for v in frames["date_obs"]]
        raw = Time(iso, format="isot", scale="utc").mjd
    t, suffix = _style.relative_time(raw)
    return t, f"MJD{suffix}"


def _by_group(frames, group="unit"):
    """Yield ``(index, label, mask)`` per distinct value of ``group``."""
    if group not in frames.colnames:
        yield 0, "all", np.ones(len(frames), dtype=bool)
        return
    keys = np.asarray(frames[group]).astype(str)
    for i, key in enumerate(sorted(set(keys.tolist()))):
        yield i, (key or "?"), keys == key


def zeropoint_vs_time(frames, group="unit", ax=None):
    """Per-frame zero point against time, one series per unit.

    Tracks transparency: a smooth downward drift is haze or rising airmass, and sudden
    dips are cloud. Each physical Seestar has its own zero point, so the series are
    expected to be offset from each other -- what matters is the shape within a series,
    not the spacing between them.
    """
    fig, ax = new_axes(ax, figsize=(6.4, 3.6))
    t, label = _time_axis(frames)
    zp = np.asarray(frames["zeropoint"], dtype=float)
    for i, name, mask in _by_group(frames, group):
        ax.plot(t[mask], zp[mask], ls="none", markersize=4.5, alpha=0.85,
                **_style.cat_kw(i, name))
    ax.set_xlabel(label)
    ax.set_ylabel("zero point (mag)")
    ax.set_title("Zero point over time")
    _style.legend(ax, title=group, outside=True)
    return ax


def metric_histogram(frames, column="rms", threshold=None, ax=None, log=False):
    """Distribution of a per-frame quality metric, with its threshold marked.

    For ``rms`` the marked line is the photometric-grade cut (0.06 mag) and the panel
    answers "how much of this dataset is usable". For ``chi2_red`` expect values of
    100-200 even for excellent frames: that is a real ~0.03 mag systematic floor rather
    than a broken fit, so read the *shape* and pick out the frames far to the right,
    rather than comparing against 1.
    """
    fig, ax = new_axes(ax, figsize=(4.6, 3.4))
    values = np.asarray(frames[column], dtype=float)
    values = values[np.isfinite(values)]
    bins = np.geomspace(max(values.min(), 1e-3), values.max(), 40) if log and len(values) \
        else 40
    ax.hist(values, bins=bins, color=_style.CATEGORICAL[0],
            edgecolor=_style.SURFACE, lw=0.6)
    if log:
        ax.set_xscale("log")
    if threshold is not None:
        _style.reference_line(ax, threshold, axis="x", label=f"{threshold:g}")
        _style.annotate(
            ax,
            f"below {threshold:g}: {int((values < threshold).sum())}/{len(values)}",
            loc="upper right",
        )
    ax.set_xlabel(column)
    ax.set_ylabel("frames")
    ax.set_title(f"Distribution of {column}")
    return ax


def depth_vs_exptime(frames, group="model", ax=None):
    """5-sigma limit against on-sky integration, with the sqrt(t) expectation drawn.

    Background-limited photometry deepens by 1.25 mag per dex of exposure, so the points
    should follow that slope. Flattening at long exposures means a systematic floor has
    taken over and the theoretical scaling no longer applies -- the caveat that stops
    the sqrt(t) extrapolation being quoted for multi-hour stacks.
    """
    fig, ax = new_axes(ax)
    t = np.asarray(frames["total_exptime"], dtype=float)
    v = np.asarray(frames["v_lim_5sigma"], dtype=float)
    ok = np.isfinite(t) & (t > 0) & np.isfinite(v)
    for i, name, mask in _by_group(frames, group):
        m = mask & ok
        ax.plot(t[m], v[m], ls="none", markersize=4.5, alpha=0.8,
                **_style.cat_kw(i, name))
    if ok.sum() > 2:
        t_ref = float(np.median(t[ok]))
        v_ref = float(np.median(v[ok]))
        grid = np.geomspace(t[ok].min(), t[ok].max(), 20)
        ax.plot(grid, v_ref + 1.25 * np.log10(grid / t_ref),
                color=_style.REFERENCE, ls="--", lw=1.1,
                label=r"$\sqrt{t}$ (1.25 mag/dex)")
    ax.set_xscale("log")
    ax.set_xlabel("total on-sky exposure (s)")
    ax.set_ylabel(r"5$\sigma$ limiting V (mag)")
    ax.set_title("Depth vs integration time")
    _style.mag_axis(ax, "y")
    _style.legend(ax, title=group, outside=True)
    return ax


def depth_vs_driver(frames, driver="sky_sb", ax=None):
    """5-sigma limit against a condition that drives it (sky brightness, seeing).

    Depth is set by integration *and* conditions. Exposure is divided out first by
    rescaling every frame to a common 15 minutes, so what is left is the conditions
    trend -- darker sky and tighter PSF both go deeper. This is the relation the
    per-unit depth model fits, and it is worth checking the frames actually span a
    range of the driver before believing a fitted coefficient.
    """
    from .depth import EXPTIME_REF, scale_limit

    fig, ax = new_axes(ax)
    x = np.asarray(frames[driver], dtype=float)
    v = scale_limit(
        np.asarray(frames["v_lim_5sigma"], dtype=float),
        np.asarray(frames["total_exptime"], dtype=float), EXPTIME_REF,
    )
    ok = np.isfinite(x) & np.isfinite(v)
    for i, name, mask in _by_group(frames, "unit"):
        m = mask & ok
        ax.plot(x[m], v[m], ls="none", markersize=4.5, alpha=0.8,
                **_style.cat_kw(i, name))
    if ok.sum() > 5:
        slope, intercept = np.polyfit(x[ok], v[ok], 1)
        grid = np.linspace(x[ok].min(), x[ok].max(), 10)
        ax.plot(grid, slope * grid + intercept, color=_style.MODEL, lw=1.3,
                label=f"slope {slope:+.2f}")
    labels = {"sky_sb": "sky brightness (mag/arcsec$^2$)",
              "fwhm_G": "green FWHM (pixels)",
              "airmass": "airmass"}
    ax.set_xlabel(labels.get(driver, driver))
    ax.set_ylabel(r"5$\sigma$ limiting V at 15 min (mag)")
    ax.set_title(f"Depth vs {driver}")
    _style.mag_axis(ax, "y")
    _style.legend(ax, outside=True)
    return ax


def condition_vs_time(frames, column="fwhm_G", ax=None):
    """An observing condition across the night, one series per unit.

    Seeing, airmass and sky brightness all evolve through a session; seeing this makes
    a feature in a light curve interpretable (or dismissable) rather than mysterious.
    """
    fig, ax = new_axes(ax, figsize=(6.4, 3.0))
    t, label = _time_axis(frames)
    y = np.asarray(frames[column], dtype=float)
    for i, name, mask in _by_group(frames, "unit"):
        ax.plot(t[mask], y[mask], ls="none", markersize=4.0, alpha=0.8,
                **_style.cat_kw(i, name))
    labels = {"fwhm_G": "green FWHM (pixels)", "airmass": "airmass",
              "sky_sb": "sky brightness (mag/arcsec$^2$)",
              "ccd_temp": "sensor temperature (C)"}
    ax.set_xlabel(label)
    ax.set_ylabel(labels.get(column, column))
    ax.set_title(f"{column} over time")
    if column == "sky_sb":
        _style.mag_axis(ax, "y")
    _style.legend(ax, title="unit", outside=True)
    return ax


def calibration_coverage(frames, ax=None):
    """Number of calibration stars per frame, over time.

    A frame whose ``n_cal`` collapses has a zero point fit from a handful of stars and
    should be distrusted whatever its rms says -- with few stars the sigma clipping has
    nothing to work with, so the scatter can look deceptively small.
    """
    fig, ax = new_axes(ax, figsize=(6.4, 3.0))
    t, label = _time_axis(frames)
    n_cal = np.asarray(frames["n_cal"], dtype=float)
    for i, name, mask in _by_group(frames, "unit"):
        ax.plot(t[mask], n_cal[mask], ls="none", markersize=4.0, alpha=0.8,
                **_style.cat_kw(i, name))
    _style.reference_line(ax, 20, axis="y", label="sparse")
    ax.set_xlabel(label)
    ax.set_ylabel("calibration stars")
    ax.set_title("Calibration coverage per frame")
    _style.legend(ax, title="unit", outside=True)
    return ax


# --- light curve ----------------------------------------------------------------------

def finder_chart(frame, wcs, stars, comps=None, target_id=None, ax=None):
    """The field with the target and the chosen comparisons marked.

    Confirms the ensemble is what you meant to select: that the target really is the
    target, and that the comparisons are spread around it rather than clustered on one
    side (which would leave a spatial systematic uncorrected).
    """
    fig, ax = new_axes(ax, figsize=(5.6, 5.6))
    green = np.asarray(frame.g, dtype=float)
    lo, hi = np.nanpercentile(green, [30, 99.5])
    ax.imshow(green, origin="lower", cmap="gray_r", vmin=lo, vmax=hi)

    def _xy(rows):
        return wcs.world_to_pixel_values(
            np.asarray(rows["ra"], dtype=float), np.asarray(rows["dec"], dtype=float)
        )

    if comps is not None and len(comps):
        cx, cy = _xy(comps)
        ax.scatter(cx, cy, s=150, facecolors="none",
                   edgecolors=_style.CATEGORICAL[0], linewidths=1.3,
                   label=f"comparisons ({len(comps)})")
        mag_col = "v_jkc_mag" if "v_jkc_mag" in comps.colnames else None
        if mag_col:
            for xi, yi, m in zip(cx, cy, np.asarray(comps[mag_col], dtype=float)):
                ax.annotate(f"{m:.1f}", (xi, yi), textcoords="offset points",
                            xytext=(9, 5), fontsize=6.5, color=_style.INK_SECONDARY)

    target = stars[np.asarray(stars["source_id"]) == target_id] if target_id is not None \
        else stars[np.asarray(stars["is_target"])]
    if len(target):
        tx, ty = _xy(target)
        ax.scatter(tx, ty, s=220, facecolors="none",
                   edgecolors=_style.CATEGORICAL[1], linewidths=1.8, label="target")
    ax.set_title("Finder chart")
    ax.set_xlabel("x (pixels)")
    ax.set_ylabel("y (pixels)")
    ax.grid(visible=False)
    ax.legend(loc="upper right")
    return ax


def lightcurve(lc, mag_col="dmag", ax=None, title=None):
    """The differential light curve: magnitude against time, with error bars.

    The headline result. The quoted scatter is the honest precision figure -- compare it
    against the median error bar: scatter much larger than the errors means a
    systematic remains, while agreement means you are at the photon limit.
    """
    fig, ax = new_axes(ax, figsize=(6.6, 3.6))
    t, suffix = _style.relative_time(np.asarray(lc["time"], dtype=float))
    y = np.asarray(lc[mag_col], dtype=float)
    err = np.asarray(lc["mag_err"], dtype=float) if "mag_err" in lc.colnames else None
    ax.errorbar(t, y, yerr=err, ls="none", marker="o", markersize=4.0,
                color=_style.BAND_COLOR.get(lc.meta.get("band", "G"), _style.CATEGORICAL[0]),
                ecolor=_style.AXIS, elinewidth=0.8, capsize=0, alpha=0.9)
    _style.reference_line(ax, 0.0 if mag_col == "dmag" else float(np.nanmedian(y)),
                          axis="y")
    ax.set_xlabel(f"BJD (TDB){suffix}")
    ax.set_ylabel("$\\Delta$V (mag)" if mag_col == "dmag" else "V (mag)")
    ax.set_title(title or f"Differential light curve, {lc.meta.get('band', '?')} band")
    _style.mag_axis(ax, "y")
    scatter = float(np.nanstd(y))
    lines = [f"n = {len(lc)}", f"scatter = {scatter * 1000:.1f} mmag"]
    if err is not None and np.isfinite(err).any():
        lines.append(f"median error = {np.nanmedian(err) * 1000:.1f} mmag")
    _style.annotate(ax, "\n".join(lines))
    return ax


def raw_flux(measurements, source_id, band="G", ax=None):
    """Uncalibrated target flux against time, beside what the ensemble removes.

    Raw flux swings with transparency and airmass by far more than any real signal, so
    a raw curve that looks like noise and a differential curve that looks clean is the
    expected outcome -- and seeing the raw version is how you confirm the differential
    step is doing real work rather than flattening a genuine signal.
    """
    fig, ax = new_axes(ax, figsize=(6.6, 3.0))
    m = measurements[
        (np.asarray(measurements["band"]) == band)
        & (np.asarray(measurements["source_id"]) == source_id)
    ]
    t, suffix = _style.relative_time(np.asarray(m["bjd_tdb"], dtype=float))
    flux = np.asarray(m["flux"], dtype=float)
    order = np.argsort(t)
    ax.plot(t[order], flux[order], ls="none", marker="o", markersize=3.5,
            color=_style.INK_MUTED)
    ax.set_xlabel(f"BJD (TDB){suffix}")
    ax.set_ylabel("aperture flux (ADU)")
    ax.set_title("Raw target flux (before differential correction)")
    return ax


def ensemble_zeropoint(lc, ax=None):
    """The per-frame ensemble zero point, with the spread across comparisons.

    The correction being applied, made visible. The band is the star-to-star scatter of
    the per-comparison zero points; it should stay roughly constant. A frame where it
    balloons had a comparison misbehave, and the target's point in that frame inherits
    the problem.
    """
    fig, ax = new_axes(ax, figsize=(6.6, 3.0))
    t, suffix = _style.relative_time(np.asarray(lc["time"], dtype=float))
    zp = np.asarray(lc["zp_ens"], dtype=float)
    spread = np.asarray(lc["zp_scatter"], dtype=float)
    ax.fill_between(t, zp - spread, zp + spread, color=_style.CATEGORICAL[0],
                    alpha=0.16, lw=0, label="comparison scatter")
    ax.plot(t, zp, ls="none", marker="o", markersize=3.5,
            color=_style.CATEGORICAL[0], label="ensemble ZP")
    ax.set_xlabel(f"BJD (TDB){suffix}")
    ax.set_ylabel("ensemble zero point (mag)")
    ax.set_title("Per-frame ensemble zero point")
    _style.legend(ax)
    return ax


def comparison_grid(curves, comps=None, ncols=4, max_panels=16):
    """A small-multiples grid: each comparison's own differential curve.

    The single most informative check on an ensemble, and worth reading before the
    target's curve. Every panel should be flat at the noise floor. A panel with a
    trend, a step or obvious excess scatter is a comparison that does not belong --
    left in, it inflates the target's error bars and can imprint its own variability
    on the target.

    Panels are ordered worst-scatter-first, so problems are in the top-left. The
    y-scale is shared so the panels are directly comparable.
    """
    import matplotlib.pyplot as plt

    ranked = sorted(curves.items(), key=lambda kv: -_scatter_of(kv[1]))
    shown = ranked[:max_panels]
    n = len(shown)
    if not n:
        raise ValueError("no comparison curves to plot")
    nrows = int(np.ceil(n / ncols))
    with plt.rc_context(_style.rc()):
        fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.1 * nrows),
                                 sharey=True, squeeze=False)
    span = max(_scatter_of(c) for _, c in shown) * 4 or 0.05
    # One epoch across every panel, so panels are directly comparable and the x axis
    # carries a readable number instead of a shared offset annotation.
    all_t = np.concatenate([np.asarray(c["time"], dtype=float)
                            for _, c in shown if len(c)])
    _, suffix = _style.relative_time(all_t)
    epoch = float(np.floor(np.nanmin(all_t)))
    for idx, (ax, (sid, lc)) in enumerate(zip(axes.ravel(), shown)):
        if not len(lc):
            ax.set_axis_off()
            continue
        scatter = _scatter_of(lc)
        # Flag the worst offenders with a reserved status colour *and* the word CHECK,
        # so the state is never carried by hue alone.
        bad = scatter > 0.05
        ax.plot(np.asarray(lc["time"], dtype=float) - epoch,
                np.asarray(lc["dmag"], dtype=float),
                ls="none", marker="o", markersize=2.6,
                color=_style.STATUS["critical"] if bad else _style.CATEGORICAL[0],
                alpha=0.85)
        ax.axhline(0.0, color=_style.AXIS, lw=0.8)
        ax.set_ylim(-span, span)
        mag = _comp_mag(comps, sid)
        title = f"{sid}" if mag is None else f"V={mag:.1f}"
        ax.set_title(f"{title}   {scatter * 1000:.0f} mmag"
                     + ("  CHECK" if bad else ""), fontsize=8)
        ax.tick_params(labelsize=6.5)
        # Tick labels on the outer edges only: inside a grid they collide with the
        # title of the panel below.
        if idx // ncols < nrows - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel(f"BJD{suffix}", fontsize=7.5)
        if idx % ncols == 0:
            ax.set_ylabel("$\\Delta$mag", fontsize=7.5)
    for ax in axes.ravel()[n:]:
        ax.set_axis_off()
    fig.suptitle("Each comparison star measured against the others", fontsize=10)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    if len(ranked) > max_panels:
        fig.text(0.5, 0.005,
                 f"showing the {max_panels} worst of {len(ranked)} comparisons",
                 ha="center", fontsize=7.5, color=_style.INK_SECONDARY)
    return fig, axes


def _scatter_of(lc):
    if not len(lc):
        return 0.0
    value = float(np.nanstd(np.asarray(lc["dmag"], dtype=float)))
    return value if np.isfinite(value) else 0.0


def _comp_mag(comps, source_id):
    if comps is None or "v_jkc_mag" not in comps.colnames:
        return None
    hit = comps[np.asarray(comps["source_id"]) == source_id]
    return float(hit["v_jkc_mag"][0]) if len(hit) else None


def scatter_vs_magnitude(curves, comps, ax=None):
    """Achieved scatter of each comparison against its brightness.

    The noise-floor panel. Bright stars should sit on a floor set by systematics and
    faint ones should climb along the photon-noise curve. Where the floor sits *is* the
    precision of the dataset; a bright star well above the floor is a bad comparison,
    and one at the very bright end climbing again is saturating.
    """
    fig, ax = new_axes(ax)
    mags, scatters = [], []
    for sid, lc in curves.items():
        mag = _comp_mag(comps, sid)
        if mag is None or not len(lc):
            continue
        mags.append(mag)
        scatters.append(_scatter_of(lc) * 1000)
    if not mags:
        ax.text(0.5, 0.5, "no comparison magnitudes available", ha="center",
                va="center", transform=ax.transAxes, color=_style.INK_SECONDARY)
        ax.set_axis_off()
        return ax
    mags = np.asarray(mags)
    scatters = np.asarray(scatters)
    ax.plot(mags, scatters, ls="none", marker="o", markersize=6,
            color=_style.CATEGORICAL[0], label="comparison stars")
    floor = float(np.nanmin(scatters[scatters > 0])) if (scatters > 0).any() else np.nan
    if np.isfinite(floor):
        _style.reference_line(ax, floor, axis="y", label=f"floor {floor:.0f} mmag")
    ax.set_yscale("log")
    ax.set_xlabel("catalogue V (mag)")
    ax.set_ylabel("achieved scatter (mmag)")
    ax.set_title("Precision vs comparison brightness")
    _style.legend(ax)
    return ax


def periodogram(pg, ax=None):
    """Lomb-Scargle power against period, with the best peak marked.

    Check the peak against the sampling before believing it: nightly cadence puts
    strong power at 1 day and its harmonics, so a period near an integer fraction of a
    day needs the phase fold to confirm it is real.
    """
    fig, ax = new_axes(ax)
    period = 1.0 / np.asarray(pg["frequency"], dtype=float)
    ax.plot(period, np.asarray(pg["power"], dtype=float),
            color=_style.CATEGORICAL[0], lw=1.1)
    ax.axvline(pg["best_period"], color=_style.MODEL, lw=1.2, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("period (days)")
    ax.set_ylabel("Lomb-Scargle power")
    ax.set_title("Periodogram")
    _style.annotate(
        ax,
        f"best period = {pg['best_period']:.5f} d\n"
        f"power = {pg['best_power']:.3f}\n"
        f"FAP = {pg['fap']:.2e}",
        loc="upper right",
    )
    return ax


def phase_fold(lc, period, t0=None, n_bins=25, ax=None):
    """The light curve folded on a period, with a binned mean over the top.

    The confirmation panel. A real periodic signal folds into a coherent, repeatable
    shape; an alias or a spurious peak folds into scatter. Two cycles are drawn so the
    shape reads continuously across the wrap.
    """
    from .lightcurves import phase_fold as _fold

    fig, ax = new_axes(ax, figsize=(5.6, 3.6))
    phase, mag = _fold(lc, period, t0=t0)
    for offset in (0.0, 1.0):
        ax.plot(phase + offset, mag, ls="none", marker="o", markersize=3.0,
                color=_style.INK_MUTED, alpha=0.55,
                label="measurements" if offset == 0 else None)
    edges = np.linspace(0, 1, n_bins + 1)
    idx = np.digitize(phase, edges) - 1
    centres, means = [], []
    for i in range(n_bins):
        sel = idx == i
        if sel.sum():
            centres.append(0.5 * (edges[i] + edges[i + 1]))
            means.append(float(np.nanmean(mag[sel])))
    for offset in (0.0, 1.0):
        ax.plot(np.asarray(centres) + offset, means, marker="o", markersize=5,
                color=_style.CATEGORICAL[1], lw=1.6,
                label="binned mean" if offset == 0 else None)
    ax.set_xlabel("phase")
    ax.set_ylabel("$\\Delta$V (mag)")
    ax.set_title(f"Phase folded on {period:.5f} d")
    ax.set_xlim(0, 2)
    _style.mag_axis(ax, "y")
    _style.legend(ax)
    return ax


def host_cutout(cut, contamination=None, ax=None):
    """The target on its host, with the aperture and the azimuthal sampling ring.

    Shows exactly what the contamination estimate measured: the ring is the
    galactocentric radius the host flux was sampled at, and each marker one clean
    azimuth. Few markers, or markers all on one side, means the estimate rests on a
    poor sample and its scatter should be treated as a real uncertainty on the result.
    """
    from matplotlib.patches import Circle

    fig, ax = new_axes(ax, figsize=(5.0, 5.0))
    stamp = np.asarray(cut["stamp"], dtype=float)
    lo, hi = np.nanpercentile(stamp, [20, 99.7])
    ax.imshow(stamp, origin="lower", cmap="gray_r", vmin=lo, vmax=hi)
    tx, ty = cut["target_xy"]
    nx_, ny_ = cut["nucleus_xy"]
    radius = cut["radius_px"]
    ax.add_patch(Circle((nx_, ny_), radius, facecolor="none",
                        edgecolor=_style.CATEGORICAL[0], lw=1.0, ls="--"))
    ax.add_patch(Circle((tx, ty), cut["aperture"], facecolor="none",
                        edgecolor=_style.CATEGORICAL[1], lw=1.6))
    ax.plot(nx_, ny_, marker="+", markersize=10, color=_style.CATEGORICAL[0],
            label="host nucleus")
    ax.plot([], [], marker="o", ls="none", color=_style.CATEGORICAL[1],
            label="target aperture")
    ax.plot([], [], ls="--", color=_style.CATEGORICAL[0], label="sampling radius")
    ax.set_title("Target on its host")
    ax.grid(visible=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(loc="upper right")
    if contamination:
        _style.annotate(
            ax,
            f"host flux {contamination['adu']:.0f} ADU\n"
            f"scatter {contamination['std']:.0f}\n"
            f"{contamination['n_azimuth']} clean azimuths",
        )
    return ax
