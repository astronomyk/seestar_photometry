"""Shared figure styling, so every diagnostic reads as one system.

Matplotlib only, and **committed to a light surface**: these figures end up in
notebooks, papers and PNGs on disk, all of which are light. There is no dark variant
to keep in sync.

Band colours
------------
R/G/B are drawn red/green/blue because in this domain the colour *is* the band name --
relabelling them would be actively misleading to the reader. That choice is not free:
red vs green is the classic colour-vision-deficiency collision, and measuring the
chosen steps (``#e34948`` / ``#008300`` / ``#2a78d6``) gives a worst-pair CVD
separation of ΔE 7.2 (protan), inside the 6-8 floor band rather than clear of it.
Alternatives were measured too -- every greener or lighter green traded protan
separation for a worse tritan one, so these steps are the best available under the
constraint.

The floor band is only legal **with secondary encoding**, so hue is never the sole
carrier of band identity here: each band also gets its own marker shape and line
style (:data:`BAND_MARKER`, :data:`BAND_STYLE`), and multi-band panels are always
legended. Do not remove those when editing a figure -- they are the accessibility
mechanism, not decoration.

Categorical series (telescope units, datasets) use a validated palette in fixed slot
order, plus the same marker-shape channel, since a scatter with more than three
series cannot clear the all-pairs floors on hue alone.
"""

#: Band hues. Domain-mandated; see the module docstring for the CVD measurement.
BAND_COLOR = {"R": "#e34948", "G": "#008300", "B": "#2a78d6"}

#: Per-band marker shape -- the secondary encoding that makes the band trio legal.
BAND_MARKER = {"R": "o", "G": "s", "B": "^"}

#: Per-band line style -- secondary encoding for line panels.
BAND_STYLE = {"R": "-", "G": "--", "B": ":"}

#: Categorical slots for non-band series (units, datasets), in fixed order. Assigned
#: by position and never cycled: a filter that drops a series must not repaint the
#: survivors. Past the eighth, fold into "other" or facet.
CATEGORICAL = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#4a3aa7", "#e34948",
)

#: Marker shapes paired with the categorical slots, so identity survives CVD and
#: greyscale printing.
CAT_MARKER = ("o", "s", "^", "D", "v", "P", "X", "*")

# Chrome and ink. Grid and axes stay recessive so the data carries the figure.
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"

#: Neutral for reference lines, models and annotations that aren't a data series.
REFERENCE = "#52514e"

#: Emphasis for a fitted model drawn over data.
MODEL = "#4a3aa7"

#: Status colours, reserved. Never reused as a series colour, and always paired with
#: a text label so state is not carried by hue alone.
STATUS = {"good": "#0ca30c", "warning": "#fab219",
          "serious": "#ec835a", "critical": "#d03b3b"}


def band_kw(band, filled=True):
    """Colour + marker keywords for one band, with the secondary encoding applied."""
    return {
        "color": BAND_COLOR[band],
        "marker": BAND_MARKER[band],
        "label": band,
        "markerfacecolor": BAND_COLOR[band] if filled else "none",
    }


def cat_kw(index, label=None):
    """Colour + marker keywords for categorical slot ``index`` (wraps by position)."""
    i = int(index) % len(CATEGORICAL)
    kw = {"color": CATEGORICAL[i], "marker": CAT_MARKER[i]}
    if label is not None:
        kw["label"] = str(label)
    return kw


def rc():
    """Matplotlib rcParams for the diagnostic figures."""
    return {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.labelcolor": INK_SECONDARY,
        "axes.edgecolor": AXIS,
        "axes.titlecolor": INK,
        "axes.linewidth": 0.8,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "lines.linewidth": 1.6,
        "figure.dpi": 110,
        "savefig.dpi": 130,
        "savefig.bbox": "tight",
    }


def use_style():
    """Apply the style globally. Called by :mod:`report`; harmless to call twice."""
    import matplotlib.pyplot as plt

    plt.rcParams.update(rc())


def new_axes(ax, figsize=(5.2, 3.6)):
    """Return ``(fig, ax)``, creating a styled figure if ``ax`` is None."""
    import matplotlib.pyplot as plt

    if ax is not None:
        return ax.figure, ax
    with plt.rc_context(rc()):
        fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


def annotate(ax, text, loc="upper left"):
    """Put a small stats block in a corner, in muted ink rather than a series colour."""
    x, ha = (0.02, "left") if "left" in loc else (0.98, "right")
    y, va = (0.97, "top") if "upper" in loc else (0.03, "bottom")
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va,
            fontsize=7.5, color=INK_SECONDARY, linespacing=1.35)


def reference_line(ax, value, axis="y", label=None, style="--"):
    """Draw a labelled reference line (a threshold, a zero, an expected value).

    Labels are placed just *outside* the axes rather than inside. Inside, they
    reliably collide with either the data or the corner stats block -- and a
    diagnostic figure is worth least exactly when it is cluttered.
    """
    fn = ax.axhline if axis == "y" else ax.axvline
    fn(value, color=REFERENCE, lw=0.9, ls=style, zorder=1)
    if not label:
        return
    if axis == "y":
        # Just outside the right edge: no data there, and it can't reach the title.
        ax.text(1.005, value, f"{label}", transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=7, color=REFERENCE)
    else:
        # Inside the top of the axes, rotated along the line. Outside the axes it
        # would land on the title, which is worse -- panel titles carry more.
        ax.text(value, 0.98, f"{label} ", transform=ax.get_xaxis_transform(),
                ha="right", va="top", fontsize=7, color=REFERENCE, rotation=90)


def legend(ax, title=None, outside=False, ncol=None, **kwargs):
    """Place a legend that doesn't sit on the data.

    ``outside=True`` puts it in a single row above the axes -- the only placement that
    is reliably collision-free for a scatter with several series, which is most of the
    per-unit panels. Marker scale is bumped because the dense clouds are drawn with
    small translucent points that are invisible at legend size.

    Call this *after* :meth:`set_title`: an outside legend occupies the strip the title
    normally sits in, so the title is re-applied with extra padding to clear it.
    """
    handles, labels = ax.get_legend_handles_labels()
    if not labels:
        return None
    if not outside:
        return ax.legend(title=title, ncol=ncol or 1, markerscale=2.5, **kwargs)

    # The legend title is dropped deliberately: outside-above it needs a whole extra
    # line, which then collides with the panel title. Unit serials and model names are
    # self-describing, and the panel title already supplies the context.
    result = ax.legend(
        loc="lower left", bbox_to_anchor=(0.0, 1.01, 1.0, 0.14),
        mode="expand", borderaxespad=0.0, ncol=ncol or min(len(labels), 4),
        markerscale=2.0, fontsize=7.5, **kwargs,
    )
    existing = ax.get_title()
    if existing:
        ax.set_title(existing, pad=18)
    return result


def relative_time(t):
    """``(t - epoch, axis_label_suffix)`` for a time axis.

    Matplotlib's default offset notation (a bare ``+2.46e6`` tucked under the axis) is
    easy to miss and easy to misread. Subtracting a round epoch and naming it in the
    axis label is unambiguous.
    """
    import numpy as np

    t = np.asarray(t, dtype=float)
    finite = t[np.isfinite(t)]
    if not len(finite):
        return t, ""
    epoch = float(np.floor(finite.min()))
    return t - epoch, f" - {epoch:.0f}"


def mag_axis(ax, axis="y"):
    """Invert a magnitude axis, since brighter is a smaller number."""
    (ax.invert_yaxis if axis == "y" else ax.invert_xaxis)()
