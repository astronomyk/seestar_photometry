#!/usr/bin/env python
"""Recover R/G/B quantum-efficiency curves from a sensor datasheet plot.

One digitiser for all three Seestar sensors, because the alternative is three
hand-tuned copies that drift apart. Each plot is described by a :class:`Plot` entry
saying how to find its axes and what they mean; everything after that is shared.

The curves are the only saturated colours in these figures, so a colour mask per
channel recovers each one to about a pixel. Two things make that harder than it sounds
and both are handled here:

* **Axis calibration differs per figure.** Sony's IMX585 plot has a full grid, which is
  self-calibrating: find the gridlines and the outermost ones are the axis limits. The
  OmniVision plot has had its grid cleaned off, leaving only the axis lines and their
  tick marks, so it is calibrated from the ticks instead.
* **Curves occlude each other.** Where all three converge in the NIR only the
  last-drawn one is visible; where a screenshot's mouse cursor sits on a curve there is
  simply a hole. Both come back as gaps, and gaps are interpolated across rather than
  guessed at -- with a loud failure if a channel is missing more than
  :data:`MAX_GAP_FRACTION` of its range.

Usage::

    uv run python tools/digitise_sensor_qe.py            # writes the ECSV + a check figure

Maintenance script, not shipped in the wheel. Its output feeds
``tools/build_unit_passbands.py``; see ``docs/background-info/`` for what the curves
are then used for.
"""

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

HERE = Path(__file__).resolve().parent
IMAGES = HERE / "sensor_qe"
OUT_ECSV = HERE / "sensor_qe.ecsv"
OUT_FIG = HERE / "sensor_qe" / "digitisation_check.png"

#: Output sampling, in nm. Wider than any passband needs so the IR-cut model has room.
WL_GRID = np.arange(330.0, 1021.0, 2.0)

#: Refuse a channel that is invisible over more than this fraction of the plot's width.
#: Occlusion and cursors cost a few percent; a broken colour mask costs most of it, and
#: the difference matters because both look like "some NaNs" from the inside.
MAX_GAP_FRACTION = 0.25


@dataclass
class Plot:
    """How to read one datasheet figure."""

    name: str
    image: str
    #: Axis ranges the calibration maps onto, in nm and in QE percent.
    wl_range: tuple
    qe_range: tuple
    #: "grid" -- outermost detected gridlines are the axis limits.
    #: "ticks" -- the axis lines plus their tick marks give the scale.
    calibration: str
    #: Expected (n_columns, n_rows) for "grid"; the count is the check that it worked.
    grid_shape: tuple = ()
    #: Fraction of a row or column that must be grid-coloured to count as a gridline.
    #: Per figure because a soft JPEG grid never reaches the fraction a crisp PNG does,
    #: and the count check above is what says the chosen value was right.
    grid_fraction: float = 0.75
    #: For "ticks": the nm step between x ticks, and the percent step between y ticks.
    tick_step: tuple = ()
    #: Per-channel RGB thresholds, as (channel_high, other_1_low, other_2_low).
    masks: dict = field(default_factory=dict)
    #: Channels that vanish under another in the NIR, mapped to the channel to read
    #: instead. Stated per plot rather than inferred: which curve ends up on top is a
    #: property of the draw order in that particular figure, and getting it wrong
    #: substitutes a curve that is *not* coincident.
    nir_fill: dict = field(default_factory=dict)
    nir_from: float = 900.0
    notes: str = ""


#: Sony's own masks are looser because the JPEG is soft; the two PNGs are crisp.
_SONY_MASKS = {"R": (150, 110, 110), "G": (130, 130, 130), "B": (140, 110, 130)}
_CRISP_MASKS = {"R": (120, 100, 100), "G": (110, 130, 130), "B": (120, 110, 130)}

PLOTS = (
    Plot(
        name="imx462", image="imx462_qe.jpg",
        wl_range=(400.0, 1000.0), qe_range=(0.0, 100.0),
        calibration="grid", grid_shape=(13, 11), grid_fraction=0.55,
        masks=_SONY_MASKS,
        nir_fill={"R": "B"},
        notes="Sony IMX462, the S50. Digitised in photometry-adventures-2 first; kept "
              "here so all three units come from one code path.",
    ),
    Plot(
        name="imx585", image="imx585_qe.png",
        wl_range=(400.0, 1000.0), qe_range=(0.0, 100.0),
        calibration="grid", grid_shape=(13, 11), masks=_CRISP_MASKS,
        nir_fill={"R": "B", "G": "B"}, nir_from=830.0,
        notes="Sony IMX585, the S30 Pro, from the ZWO ASI585 datasheet plot. The "
              "screenshot caught a mouse cursor over the green peak near 520 nm; it "
              "occludes the curve rather than recolouring it, so it reads as a gap.",
    ),
    Plot(
        name="os08b10", image="os08b10_qe.png",
        wl_range=(300.0, 1050.0), qe_range=(0.0, 90.0),
        calibration="ticks", tick_step=(100.0, 10.0), masks=_CRISP_MASKS,
        nir_fill={"R": "B", "G": "B"}, nir_from=830.0,
        notes="OmniVision OS08B10, the S50 Pro. The grid was cleaned off this "
              "screenshot, so it is calibrated from the axis tick marks.",
    ),
)


def _clusters(index, gap=3):
    """Mean position of each run of consecutive indices."""
    if not len(index):
        return np.array([])
    groups = np.split(index, np.where(np.diff(index) > gap)[0] + 1)
    return np.array([group.mean() for group in groups])


def _calibrate_grid(gray, shape, fraction):
    """``(x0, x1, y0, y1)`` in pixels from the outermost gridlines."""
    height, width = gray.shape
    columns = _clusters(np.flatnonzero(gray.sum(axis=0) / height > fraction))
    rows = _clusters(np.flatnonzero(gray.sum(axis=1) / width > fraction))
    if (len(columns), len(rows)) != shape:
        raise RuntimeError(
            f"gridline detection found {len(columns)} x {len(rows)} lines, "
            f"expected {shape[0]} x {shape[1]}")
    return columns[0], columns[-1], rows[0], rows[-1]


def _calibrate_ticks(dark):
    """``(x0, x1, y0, y1)`` in pixels from the axis lines and their tick marks.

    Fitted over every tick rather than taken from the first and last: a single
    mis-detected tick then moves the scale by a fraction of a tick instead of setting it.
    """
    height, width = dark.shape
    axis_x = _clusters(np.flatnonzero(dark.sum(axis=0) / height > 0.30))
    axis_y = _clusters(np.flatnonzero(dark.sum(axis=1) / width > 0.30))
    if len(axis_x) != 1 or len(axis_y) != 1:
        raise RuntimeError(f"expected one axis line each way, got "
                           f"{len(axis_x)} vertical and {len(axis_y)} horizontal")
    x_axis, y_axis = float(axis_x[0]), float(axis_y[0])

    # Ticks sit just inside the y axis and just below the x axis.
    strip = dark[:, int(x_axis) - 8:int(x_axis)]
    y_ticks = _clusters(np.flatnonzero(strip.sum(axis=1) >= strip.shape[1] * 0.6))
    y_ticks = y_ticks[y_ticks <= y_axis + 1]
    strip = dark[int(y_axis):int(y_axis) + 9, :]
    x_ticks = _clusters(np.flatnonzero(strip.sum(axis=0) >= strip.shape[0] * 0.6))
    x_ticks = x_ticks[x_ticks >= x_axis - 1]
    if len(x_ticks) < 3 or len(y_ticks) < 3:
        raise RuntimeError(f"too few ticks: {len(x_ticks)} on x, {len(y_ticks)} on y")

    return x_axis, _tick_spacing(x_ticks), y_axis, _tick_spacing(y_ticks)


def _tick_spacing(ticks):
    """Pixels per tick, from tick positions that may have gaps in them.

    Fitting position against *index* looks like the robust choice and is the opposite:
    one undetected tick shifts every index after it, so the slope absorbs the gap. On
    the OmniVision x axis, where two of eight ticks are missed, that read 192 px per
    100 nm against a true 137 -- a 40% wavelength stretch that still produced a
    plausible-looking curve.

    So the spacing is measured between *adjacent* ticks instead. A gap shows up as a
    near-integer multiple of the common step, which is divided back out; the result is
    the total span over the total number of steps it spans, which is every tick's worth
    of leverage and none of the vulnerability.
    """
    ticks = np.sort(np.asarray(ticks, dtype=float))
    steps = np.diff(ticks)
    counts = np.maximum(np.round(steps / np.median(steps)), 1.0)
    return float(steps.sum() / counts.sum())


#: Shortest run of masked pixels that can be a curve rather than a speck. The strokes
#: in these figures are 3-5 px; two is enough to keep anti-aliased edges and enough to
#: reject the single pixels left behind where a grid was erased.
MIN_STROKE = 2


def _stroke(rows):
    """The longest run of consecutive masked rows in one column, or ``None``.

    Taking the median of *every* masked pixel in the column looks equivalent and is not.
    Where a curve is hidden under another, what survives the mask is not the curve but
    whatever else happens to be that colour -- on the OmniVision screenshot, specks left
    along the axis where the grid was erased. Those are near QE = 0, so the column reads
    zero and the curve acquires a spike straight down to the axis. Every such spike in
    that figure sat above 800 nm, where all three curves overlap.
    """
    if len(rows) < MIN_STROKE:
        return None
    runs = np.split(rows, np.flatnonzero(np.diff(rows) > 1) + 1)
    longest = max(runs, key=len)
    return longest if len(longest) >= MIN_STROKE else None


def digitise(plot, images=IMAGES):
    """``(wavelength_nm, {channel: qe_fraction})`` for one datasheet figure."""
    from PIL import Image

    img = np.asarray(Image.open(images / plot.image).convert("RGB")).astype(int)
    height, width, _ = img.shape
    red, green, blue = img[..., 0], img[..., 1], img[..., 2]
    gray = (abs(red - green) < 40) & (abs(green - blue) < 40) & (red < 215)
    dark = (red < 150) & (green < 150) & (blue < 150)

    if plot.calibration == "grid":
        x0, x1, y0, y1 = _calibrate_grid(gray, plot.grid_shape, plot.grid_fraction)
        def to_wl(col):
            return (plot.wl_range[0]
                    + (col - x0) / (x1 - x0) * (plot.wl_range[1] - plot.wl_range[0]))
        def to_qe(row):
            return (plot.qe_range[1]
                    + (row - y0) / (y1 - y0) * (plot.qe_range[0] - plot.qe_range[1]))
        top, bottom = int(min(y0, y1)), int(max(y0, y1))
        columns = np.arange(int(min(x0, x1)) + 2, int(max(x0, x1)) - 1)
    else:
        x_axis, px_per_x, y_axis, px_per_y = _calibrate_ticks(dark)
        def to_wl(col):
            return plot.wl_range[0] + (col - x_axis) / px_per_x * plot.tick_step[0]
        def to_qe(row):
            return (y_axis - row) / px_per_y * plot.tick_step[1]
        top = int(y_axis - px_per_y * (plot.qe_range[1] / plot.tick_step[1]) - 6)
        bottom = int(y_axis)
        columns = np.arange(int(x_axis) + 2,
                            int(x_axis + px_per_x * ((plot.wl_range[1]
                                                      - plot.wl_range[0])
                                                     / plot.tick_step[0])))
    columns = columns[(columns >= 0) & (columns < width)]
    top, bottom = max(top, 0), min(bottom, height - 1)

    masks = {}
    for channel, (high, lo1, lo2) in plot.masks.items():
        if channel == "R":
            masks[channel] = (red > high) & (green < lo1) & (blue < lo2)
        elif channel == "G":
            masks[channel] = (green > high) & (red < lo1) & (blue < lo2)
        else:
            masks[channel] = (blue > high) & (red < lo1) & (green < lo2)

    wl = to_wl(columns.astype(float))
    curves = {}
    for channel, mask in masks.items():
        values = np.full(len(columns), np.nan)
        for i, column in enumerate(columns):
            rows = _stroke(np.flatnonzero(mask[top:bottom + 1, column]))
            if rows is not None:
                # Median of the stroke, so a 4 px line reads as its centre.
                values[i] = to_qe(top + np.median(rows))
        curves[channel] = values

    for channel, donor in plot.nir_fill.items():
        # Drawn underneath another curve where they coincide, so there is nothing to
        # read -- but they coincide, which is what makes the substitution legitimate.
        # Without it the gap runs to the edge of the plot with no anchor beyond it, and
        # the interpolation below holds the last value: on the OmniVision figure that
        # left green flat at 33% from 965 nm out to 1050 while the other two descended
        # to 10%. Harmless for photometry, since the IR cut has killed everything by
        # 700 nm -- but visibly wrong, and a curve that is visibly wrong somewhere is
        # not evidence for itself anywhere.
        hidden = np.isnan(curves[channel]) & (wl > plot.nir_from)
        curves[channel][hidden] = curves[donor][hidden]

    for channel, values in curves.items():
        missing = np.isnan(values)
        if missing.mean() > MAX_GAP_FRACTION:
            raise RuntimeError(
                f"{plot.name} {channel}: {100*missing.mean():.0f}% of the plot width has "
                "no pixels of that colour -- the mask or the calibration is wrong")
        ok = ~missing
        curves[channel] = np.clip(
            np.interp(np.arange(len(values)), np.flatnonzero(ok), values[ok]), 0.0, None)
        if missing.any():
            print(f"  {plot.name} {channel}: interpolated across "
                  f"{100*missing.mean():.1f}% of the width", flush=True)
    return wl, curves


def on_grid(plot, wl, curves):
    """The curves resampled onto :data:`WL_GRID`, as a QE fraction, edges held flat.

    Held rather than extrapolated to zero outside the plotted range: below 400 nm on the
    Sony figures the sensor response is still rising and the value at the edge is the
    only honest estimate. It does not matter much either way -- the UV cut is at 390 nm
    and kills everything down there -- but a linear extrapolation of a rising curve
    would invent response the datasheet does not show.
    """
    out = {}
    for channel, values in curves.items():
        out[channel] = np.interp(WL_GRID, wl, values / 100.0,
                                 left=values[0] / 100.0, right=values[-1] / 100.0)
    return out


def overlay_figure(plot, out_dir, images=IMAGES):
    """Draw the extracted curves back onto the source plot, in its own pixel space.

    The check that actually catches things. A digitised curve plotted on clean axes
    always looks like a QE curve, so it cannot tell you whether the *axis calibration*
    was right -- the OmniVision plot was briefly read 40% too wide in wavelength and the
    result was perfectly smooth and completely wrong. Reprojected onto the original
    image, a bad scale shows up immediately as a trace that drifts off its curve.

    The detected axis and tick positions are drawn too, since that is what the mapping
    is built from.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    from seestar_photometry import _style

    img = np.asarray(Image.open(images / plot.image).convert("RGB"))
    wl, curves = digitise(plot, images)
    to_col, to_row, ticks = _pixel_mapping(plot, images)

    with plt.rc_context(_style.rc()):
        fig, axes = plt.subplots(2, 1, figsize=(9.0, 9.4),
                                 gridspec_kw={"height_ratios": [1.35, 1.0]})
        axes[0].imshow(img)
        for x in ticks["x"]:
            axes[0].axvline(x, color="0.35", lw=0.5, ls=":", alpha=0.8)
        for y in ticks["y"]:
            axes[0].axhline(y, color="0.35", lw=0.5, ls=":", alpha=0.8)
        for channel in ("R", "G", "B"):
            axes[0].plot(to_col(wl), to_row(curves[channel]), color="k", lw=1.0,
                         ls=(0, (4, 3)), label="extracted" if channel == "R" else None)
        axes[0].set_xlim(0, img.shape[1])
        axes[0].set_ylim(img.shape[0], 0)
        axes[0].set_axis_off()
        axes[0].set_title(f"{plot.name}: extraction (black dashes) over {plot.image}, "
                          "with the detected grid", loc="left")
        axes[0].legend(loc="upper right")

        for channel in ("R", "G", "B"):
            axes[1].plot(wl, curves[channel], color=_style.BAND_COLOR[channel],
                         ls=_style.BAND_STYLE[channel], lw=1.6, label=channel)
        axes[1].set_xlim(*plot.wl_range)
        axes[1].set_ylim(0, plot.qe_range[1])
        axes[1].set_xlabel("wavelength [nm]")
        axes[1].set_ylabel("QE [%]")
        axes[1].set_title("what was extracted, on its own axes", loc="left")
        _style.legend(axes[1], ncol=3)
        fig.tight_layout()

        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"{plot.name}_extraction_check.png"
        fig.savefig(path, dpi=140)
        plt.close(fig)
    return path


def _pixel_mapping(plot, images=IMAGES):
    """``(wl -> column, qe -> row, {"x": ticks, "y": ticks})`` for one figure.

    The same calibration :func:`digitise` uses, run backwards, so the overlay is drawn
    with the mapping under test rather than with a second one that could agree by luck.
    """
    from PIL import Image

    img = np.asarray(Image.open(images / plot.image).convert("RGB")).astype(int)
    red, green, blue = img[..., 0], img[..., 1], img[..., 2]
    gray = (abs(red - green) < 40) & (abs(green - blue) < 40) & (red < 215)
    dark = (red < 150) & (green < 150) & (blue < 150)

    if plot.calibration == "grid":
        height, width = gray.shape
        xs = _clusters(np.flatnonzero(gray.sum(axis=0) / height > plot.grid_fraction))
        ys = _clusters(np.flatnonzero(gray.sum(axis=1) / width > plot.grid_fraction))
        x0, x1, y0, y1 = xs[0], xs[-1], ys[0], ys[-1]
        lo, hi = plot.wl_range
        qlo, qhi = plot.qe_range
        return (lambda w: x0 + (np.asarray(w) - lo) / (hi - lo) * (x1 - x0),
                lambda q: y0 + (qhi - np.asarray(q)) / (qhi - qlo) * (y1 - y0),
                {"x": xs, "y": ys})

    x_axis, px_per_x, y_axis, px_per_y = _calibrate_ticks(dark)
    n_x = int(round((plot.wl_range[1] - plot.wl_range[0]) / plot.tick_step[0]))
    n_y = int(round(plot.qe_range[1] / plot.tick_step[1]))
    return (
        lambda w: x_axis + (np.asarray(w) - plot.wl_range[0]) / plot.tick_step[0]
        * px_per_x,
        lambda q: y_axis - np.asarray(q) / plot.tick_step[1] * px_per_y,
        {"x": [x_axis + i * px_per_x for i in range(n_x + 1)],
         "y": [y_axis - i * px_per_y for i in range(n_y + 1)]},
    )


def main():
    import argparse

    from astropy.table import Table

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--overlay", type=Path, default=None,
                    help="also write per-figure extraction-check overlays here")
    ap.add_argument("--only", default=None, help="restrict the overlays to one plot")
    args = ap.parse_args()

    table = Table({"wavelength_nm": WL_GRID})
    meta = []
    for plot in PLOTS:
        print(f"{plot.name} <- {plot.image}", flush=True)
        wl, curves = digitise(plot)
        print(f"  {plot.wl_range[0]:.0f}-{plot.wl_range[1]:.0f} nm over "
              f"{len(wl)} columns; peak QE "
              + ", ".join(f"{c} {np.nanmax(v):.0f}%" for c, v in sorted(curves.items())),
              flush=True)
        for channel, values in on_grid(plot, wl, curves).items():
            table[f"{channel}_{plot.name}"] = values
        meta.append(f"{plot.name}: {plot.notes}")

    table.meta["comments"] = [
        "Sensor quantum efficiency, digitised from datasheet plots by colour-masking,",
        "as a fraction (not percent), resampled onto a common 330-1020 nm, 2 nm grid.",
        "Values outside each figure's plotted range are held at the edge value.",
        *meta,
        f"Built by tools/{Path(__file__).name}.",
    ]
    table.write(str(OUT_ECSV), overwrite=True)
    print(f"wrote {OUT_ECSV}")
    _check_figure(table)

    if args.overlay is not None:
        for plot in PLOTS:
            if args.only and plot.name != args.only:
                continue
            print(f"wrote {overlay_figure(plot, args.overlay)}")


def _check_figure(table):
    """Overlay every digitised curve, so a bad mask is visible rather than inferred."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from seestar_photometry import _style

    wl = np.asarray(table["wavelength_nm"])
    with plt.rc_context(_style.rc()):
        fig, axes = plt.subplots(len(PLOTS), 1, figsize=(7.4, 2.4 * len(PLOTS)),
                                 sharex=True)
        for ax, plot in zip(np.atleast_1d(axes), PLOTS):
            for channel in ("R", "G", "B"):
                ax.plot(wl, 100 * np.asarray(table[f"{channel}_{plot.name}"]),
                        color=_style.BAND_COLOR[channel],
                        ls=_style.BAND_STYLE[channel], lw=1.5, label=channel)
            ax.axvspan(*plot.wl_range, color=_style.GRID, alpha=0.5, lw=0, zorder=0)
            ax.set_ylabel("QE [%]")
            ax.set_title(plot.name, loc="left")
            ax.set_ylim(0, 100)
        np.atleast_1d(axes)[0].legend(ncol=3)
        np.atleast_1d(axes)[-1].set_xlabel("wavelength [nm]")
        np.atleast_1d(axes)[-1].set_xlim(wl[0], wl[-1])
        fig.tight_layout()
        fig.savefig(OUT_FIG)
        plt.close(fig)
    print(f"wrote {OUT_FIG}")


if __name__ == "__main__":
    main()
