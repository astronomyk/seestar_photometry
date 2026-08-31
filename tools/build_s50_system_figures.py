#!/usr/bin/env python
"""Data and figures for ``docs/background-info/s50-photometric-system.md``.

Answers four questions with measurements rather than assertions:

1. where the S50's R/G/B bands sit against Johnson-Kron-Cousins and SDSS;
2. what converts an S50 magnitude into one of those systems;
3. whether that conversion needs a term beyond a linear colour term;
4. what calibrating against the S50's own bands does to the zero-point scatter of a
   real stacked frame.

Everything in 1-3 comes from integrating **the same Gaia XP spectra** through all three
sets of passbands with the same integrator, so no cross-system zero point, aperture or
standardisation enters the comparison. Question 4 uses the 20 Cep stacks the S50 system
was derived on.

Two stages, cached:

``sample``
    Integrate a spread of bulk XP files through the S50, JKC and SDSS passbands and
    join Gaia's own ``v_jkc_mag`` onto the result. Written to ``<work>/s50_system.parquet``.
``figures``
    Fit the transformations and draw the panels.

Usage::

    uv run python tools/build_s50_system_figures.py --work D:/tmp/xp_s50 \\
        --src D:/tmp/gaia_s50 --files 24

Maintenance script, not shipped in the wheel.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import build_gaia_s50 as s50  # noqa: E402
from seestar_photometry import _style, gaiadb  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
STANDARD = REPO / "tools" / "standard_passbands.ecsv"
FIGURES = REPO / "docs" / "background-info" / "figures"

#: The S50's own bands and the standard bands each is fitted against. Chosen by measured
#: effective wavelength (:func:`band_statistics`), not by name: R_S50 lands 4 nm from
#: Sloan *r* and 26 nm from Cousins *R*, so pairing it with *i* on the strength of being
#: "the red one" would be fitting across a 127 nm gap. The first entry per system is the
#: one the figures draw; ``i_SDSS`` and ``I_JKC`` are fitted only to put a number on how
#: badly a band beyond the 693 nm IR cut extrapolates.
S50_BANDS = ("R", "G", "B")
PARTNER = {"B": ("B_JKC", "g_SDSS"),
           "G": ("V_JKC", "r_SDSS"),
           "R": ("R_JKC", "r_SDSS", "i_SDSS", "I_JKC")}
STANDARD_BANDS = ("B_JKC", "V_JKC", "R_JKC", "I_JKC", "g_SDSS", "r_SDSS", "i_SDSS")

#: Colour the transformations are expressed in. An S50 user can measure this from their
#: own frames, which is the whole point -- a transformation keyed on a colour you would
#: need a different telescope to obtain is not a transformation you can apply.
COLOUR = ("B", "R")

#: Colour range the fits are quoted over: the 1st to 99.5th percentile of the sample.
#: Measured rather than chosen -- above 1.65 the sample is 0.5% of the stars and a cubic
#: fitted over a wider range does not interpolate there, it extrapolates, which is how a
#: transformation ends up confidently wrong on exactly the red stars it was meant to fix.
COLOUR_RANGE = (0.0, 1.65)

#: Polynomial degrees compared when asking whether a linear colour term is enough.
DEGREES = (1, 2, 3)


# --- stage 1: the sample ------------------------------------------------------------------

def choose_files(entries, n):
    """``n`` bulk files spread evenly through the HEALPix ordering.

    Even steps through hp8 rather than a random draw: nested HEALPix walks the sky in
    a space-filling order, so a regular stride samples every Galactic latitude instead
    of whichever ones a seed happened to pick.
    """
    step = max(len(entries) // n, 1)
    return entries[::step][:n]


def build_sample(work, src, n_files, quiet=False):
    """Integrate a spread of XP files through every system; return one table."""
    from astropy.table import Table, vstack

    entries = s50.enumerate_files(work)
    responses = s50.load_passbands()
    standard = Table.read(str(STANDARD))
    for band in STANDARD_BANDS:
        responses[band] = np.asarray(standard[band], float)

    chunks = []
    for name, lo, hi, _ in choose_files(entries, n_files):
        path = s50.local_path(name, Path(work) / "files")
        if not path.exists():
            continue
        t0 = time.time()
        ids, flux, flux_err = s50.parse_file(path)

        row = {"source_id": ids}
        for band in S50_BANDS:
            mag, _ = s50.synth_mags(flux, flux_err, s50.XP_GRID, responses[band])
            row[f"{band}_S50"] = mag
        for band in STANDARD_BANDS:
            mag, _ = s50.synth_mags(flux, flux_err, s50.XP_GRID, responses[band])
            row[band] = mag

        # Gaia's own GSPC V, so the synthesis can be checked against the scale the
        # package actually calibrates on rather than only against itself.
        gaia = _catalogue_columns(ids, src)
        row.update(gaia)
        chunks.append(Table(row))
        if not quiet:
            print(f"  {name}: {len(ids)} spectra in {time.time()-t0:.1f} s", flush=True)

    out = vstack(chunks)
    print(f"sample: {len(out)} stars from {len(chunks)} bulk files")
    return out


def _catalogue_columns(ids, src):
    """``phot_g_mean_mag`` and Gaia's ``v_jkc_mag`` for these ids, nan where absent."""
    import pyarrow.parquet as pq

    want = ("phot_g_mean_mag", "v_jkc_mag")
    out = {c: np.full(len(ids), np.nan) for c in want}
    order = np.argsort(ids, kind="stable")
    sorted_ids = ids[order]
    for pixel in np.unique(gaiadb.hpx_of(ids)):
        path = gaiadb.part_path(int(pixel), src)
        if not path.exists():
            continue
        table = pq.read_table(path, columns=["source_id", *want])
        sid = table["source_id"].to_numpy()
        where = np.clip(np.searchsorted(sorted_ids, sid), 0, len(sorted_ids) - 1)
        hit = sorted_ids[where] == sid
        target = order[where[hit]]
        for column in want:
            out[column][target] = table[column].to_numpy(zero_copy_only=False)[hit]
    return out


# --- stage 2: the fits --------------------------------------------------------------------

def clipped_polyfit(colour, delta, degree, n_sigma=3.0, passes=2):
    """Polynomial fit of ``delta`` on ``colour`` with sigma clipping.

    Returns ``(coefficients, rms, n_used, keep_mask)``. Clipped rather than plain least
    squares because the sample is real stars: a few percent are blends, emission-line
    objects or have a bad XP solution, and they would otherwise set the slope.
    """
    keep = np.isfinite(colour) & np.isfinite(delta)
    for _ in range(passes):
        coef = np.polyfit(colour[keep], delta[keep], degree)
        resid = delta - np.polyval(coef, colour)
        scatter = np.std(resid[keep])
        keep = keep & (np.abs(resid) < n_sigma * scatter)
    coef = np.polyfit(colour[keep], delta[keep], degree)
    resid = delta[keep] - np.polyval(coef, colour[keep])
    return coef, float(np.std(resid)), int(keep.sum()), keep


def transformations(sample):
    """Fit every S50 -> standard conversion at each degree. Returns a list of rows."""
    colour = np.asarray(sample[f"{COLOUR[0]}_S50"]) - np.asarray(sample[f"{COLOUR[1]}_S50"])
    inside = (colour >= COLOUR_RANGE[0]) & (colour <= COLOUR_RANGE[1])

    rows = []
    for band in S50_BANDS:
        source = np.asarray(sample[f"{band}_S50"])
        for target in PARTNER[band]:
            delta = np.asarray(sample[target]) - source
            good = inside & np.isfinite(delta) & np.isfinite(colour)
            for degree in DEGREES:
                coef, rms, n, _ = clipped_polyfit(colour[good], delta[good], degree)
                rows.append({"s50_band": band, "target": target, "degree": degree,
                             "coef": coef, "rms": rms, "n": n})
    return rows


def band_statistics():
    """Photon-weighted effective wavelength and FWHM of every band, in nm.

    The effective wavelength is ``int(S lam^2 dlam) / int(S lam dlam)`` -- the extra
    ``lam`` in both integrals is the photon weighting the magnitudes themselves use, so
    the number quoted here is the one the passband actually observes at.
    """
    from astropy.table import Table

    responses = s50.load_passbands()
    standard = Table.read(str(STANDARD))
    for band in STANDARD_BANDS:
        responses[band] = np.asarray(standard[band], float)

    wl = s50.XP_GRID
    rows = []
    for name in [f"{b}_S50" for b in S50_BANDS] + list(STANDARD_BANDS):
        curve = responses[name.replace("_S50", "") if name.endswith("_S50") else name]
        weight = curve * wl
        lam_eff = float(np.trapezoid(weight * wl, wl) / np.trapezoid(weight, wl))
        half = curve >= 0.5 * curve.max()
        rows.append({"band": name, "lam_eff": lam_eff,
                     "fwhm": float(wl[half].max() - wl[half].min())})
    return rows


def write_coefficients(out_dir, fits):
    """The fitted transformations as an ECSV, so the doc's tables are not the only copy."""
    from astropy.table import Table

    order = max(DEGREES)
    rows = {"s50_band": [], "target": [], "degree": [], "rms_mmag": [], "n": []}
    rows.update({f"c{i}": [] for i in range(order + 1)})
    for row in fits:
        rows["s50_band"].append(row["s50_band"])
        rows["target"].append(row["target"])
        rows["degree"].append(row["degree"])
        rows["rms_mmag"].append(round(1000 * row["rms"], 2))
        rows["n"].append(row["n"])
        # numpy returns highest power first; the table is written lowest first, which is
        # the order the doc writes the polynomial in and the order a reader will expect.
        ascending = list(row["coef"][::-1])
        for i in range(order + 1):
            rows[f"c{i}"].append(round(ascending[i], 6) if i < len(ascending) else None)

    table = Table(rows)
    table.meta["comments"] = [
        "S50 -> standard-system transformations, delta = target - S50_band, as",
        f"delta = sum_i c_i * x^i with x = B_S50 - R_S50 over "
        f"[{COLOUR_RANGE[0]}, {COLOUR_RANGE[1]}].",
        "All magnitudes are AB and synthesised from the same Gaia XP spectra.",
        "rms is the sigma-clipped residual scatter of the fit.",
        "Built by tools/build_s50_system_figures.py.",
    ]
    table.write(str(out_dir / "s50-transformations.ecsv"), overwrite=True)
    return table


def gaia_scale_offset(sample):
    """Offset and scatter between the synthetic AB V here and Gaia's GSPC ``v_jkc_mag``.

    The two are on different zero-point conventions -- this integrator is AB, Gaia's
    JKC is Vega -- so the mean difference is the AB-Vega offset for V and is expected.
    The *scatter* about it is the number that matters: it is how far this synthesis and
    Gaia's standardised one disagree star by star.
    """
    delta = np.asarray(sample["v_jkc_mag"]) - np.asarray(sample["V_JKC"])
    good = np.isfinite(delta)
    median = float(np.median(delta[good]))
    resid = delta[good] - median
    keep = np.abs(resid) < 3 * np.std(resid)
    return median, float(np.std(resid[keep])), int(keep.sum())


# --- stage 3: the zero point of a real frame ------------------------------------------------

#: One 15-minute 20 Cep stack, and the reference column each plane is calibrated against
#: under the two schemes. The JKC pairing is the status quo: the closest classical band.
FRAME = Path("D:/seestar_saturation_data_set/stacks_local/n2_c17_dither/X15min/"
             "n2_c17_ditherlocal_X015min_20260812-2130_n17_t340s_20s_IRCUT.fit")
REFERENCES = {"G": ("v_jkc_mag", "g_s50_mag"),
              "R": ("r_jkc_mag", "r_s50_mag"),
              "B": ("b_jkc_mag", "b_s50_mag")}

#: Calibrator cuts, matching ``photometry-adventures-2/paper/scripts/test_bayer_zp.py``
#: so the single frame here can be read against that survey's 57 stacks.
CAL_RANGE = (9.5, 13.0)
SNR_MIN = 30.0
SAT_FRAC = 0.90
CLIP_ADU = 65535.0
ENCLOSED = 0.95


def measure_frame(frame_path, src, quiet=False):
    """Forced photometry of one stack against the ``v2`` catalogue, **unfiltered**.

    Returns a per-plane dict of equal-length arrays: every catalogue source that landed
    on the chip, with the flags a calibrator cut needs rather than the cut applied. The
    zero-point-linearity figure needs the saturated and the faint ends that
    :func:`calibrator_mask` throws away, so the cut lives outside the measurement.

    The catalogue comes from :func:`gaiadb.cone` on the built ``v2`` dataset -- which is
    the point, since that dataset is what a user would have -- joined to ``bp_rp`` and
    ``phot_variable_flag`` from the field's TAP result. Those two columns are absent from
    the all-sky build (they live only in ``gaia_source``), and BP-RP in particular has to
    come from outside: a colour axis built from the S50 bands would be testing the S50
    reference against itself.
    """
    from astropy.table import Table, join
    from seestar_photometry import astrometry, catalogs, frames, photometry

    frame = frames.load_frame(str(frame_path))
    wcs = astrometry.solve(frame, solver="astap")

    ny, nx = frame.shape
    centre = wcs.pixel_to_world(nx / 2, ny / 2)
    catalogue = gaiadb.cone((centre.ra.deg, centre.dec.deg), 1.7, directory=src)
    extra = Table.read(r"D:/CrowdSky/crowdsky-workspace/Repos/photometry-adventures-2"
                       r"/paper/data/gaia_20Cep.ecsv")["source_id", "bp_rp",
                                                       "phot_variable_flag"]
    # Both columns exist in the schema but are entirely masked in the all-sky build, and
    # a join would silently rename the pair to bp_rp_1/bp_rp_2 rather than fill them.
    catalogue.remove_columns(["bp_rp", "phot_variable_flag"])
    catalogue = join(catalogue, extra, keys="source_id", join_type="inner")
    inframe = catalogs.sources_in_frame(catalogue, wcs, frame)
    if not quiet:
        print(f"  {frame_path.name}: {len(inframe)} catalogue sources on the frame")

    phot = photometry.forced_photometry(
        frame, inframe["ra"], inframe["dec"], wcs,
        source_id=inframe["source_id"], enclosed=ENCLOSED)
    blended = _blend_flag(inframe, catalogue, phot.meta["aperture"][1],
                          astrometry.pixel_scale(frame))
    if not quiet:
        print(f"  aperture {phot.meta['aperture'][1]:.2f} px (green), "
              f"{blended.sum()} of {len(inframe)} sources blended", flush=True)

    lut = {int(s): i for i, s in enumerate(inframe["source_id"])}
    out = {}
    for plane in REFERENCES:
        rows = phot[(phot["band"] == plane) & phot["on_chip"]]
        index = np.array([lut[int(s)] for s in rows["source_id"]])
        flux = np.asarray(rows["flux"], float)
        with np.errstate(invalid="ignore", divide="ignore"):
            m_inst = -2.5 * np.log10(np.where(flux > 0, flux, np.nan))

        def column(name):
            return np.ma.filled(np.asarray(inframe[name])[index], np.nan).astype(float)

        out[plane] = {
            "m_inst": m_inst,
            "bp_rp": column("bp_rp"),
            # Carried for every plane, not just G: the calibrator window is defined in
            # synthetic V whichever plane is being calibrated, so that one sample of
            # stars is shared and the planes stay comparable.
            "v_jkc_mag": column("v_jkc_mag"),
            "blend": blended[index],
            "snr": np.asarray(rows["snr"], float),
            "peak": np.asarray(rows["max_pix_value"], float),
            "variable": np.array([str(s) for s in
                                  np.asarray(inframe["phot_variable_flag"])[index]])
            == "VARIABLE",
            **{ref: column(ref) for ref in REFERENCES[plane]},
        }
    return out


#: A neighbour inside this many aperture radii and no more than ``BLEND_DMAG`` fainter
#: puts flux in the aperture that does not belong to the star. At 2.4 arcsec/px with a
#: ~10 arcsec FWHM this is the norm rather than the exception in a rich field, which is
#: why dropping the cut costs more than it looks like it should.
BLEND_APER = 1.5
BLEND_DMAG = 2.5


def _blend_flag(inframe, catalogue, aperture_px, pixel_scale):
    """True where a catalogue neighbour sits close enough to contaminate the aperture."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    radius = BLEND_APER * aperture_px * pixel_scale
    here = SkyCoord(inframe["ra"], inframe["dec"], unit="deg")
    everything = SkyCoord(catalogue["ra"], catalogue["dec"], unit="deg")
    g_all = np.asarray(catalogue["phot_g_mean_mag"], float)
    g_here = np.asarray(inframe["phot_g_mean_mag"], float)

    i_all, i_here, sep, _ = here.search_around_sky(everything, radius * u.arcsec)
    flag = np.zeros(len(inframe), bool)
    # The self-match is excluded by separation rather than by identifier: a source is its
    # own nearest neighbour at exactly zero, and nothing else lands that close.
    real = sep.arcsec > 1e-6
    contaminating = g_all[i_all[real]] < g_here[i_here[real]] + BLEND_DMAG
    flag[i_here[real][contaminating]] = True
    return flag


def calibrator_mask(data, snr_min=SNR_MIN):
    """The calibrator cut: bright enough, unsaturated, unblended, steady.

    Matches ``photometry-adventures-2/paper/scripts/test_bayer_zp.py`` so the numbers
    here can be read against that survey's 57 stacks. The blend term is the one that is
    easy to leave out and expensive to: without it this frame's calibrator rms is 164
    mmag rather than the 87 the survey reports, because a 60-minute stack of a Galactic
    field is crowded and a blended star's aperture is simply measuring two stars.
    """
    return (np.isfinite(data["m_inst"]) & np.isfinite(data["bp_rp"])
            & (data["snr"] > snr_min)
            & (data["peak"] < SAT_FRAC * CLIP_ADU)
            & (data["v_jkc_mag"] >= CAL_RANGE[0])
            & (data["v_jkc_mag"] < CAL_RANGE[1])
            & ~data["blend"] & ~data["variable"])


def frame_calibrators(frame_path, src, quiet=False):
    """:func:`measure_frame` with the calibrator cut applied, per plane."""
    out = {}
    for plane, data in measure_frame(frame_path, src, quiet).items():
        keep = calibrator_mask(data)
        out[plane] = {key: value[keep] for key, value in data.items()}
    return out


def zeropoint_fit(reference, m_inst, colour):
    """Median zero point, then the residual slope against colour.

    Deliberately *not* the package's :func:`calibration.fit_zeropoint`, which fits the
    colour term as a free parameter. The question here is how large a colour term the
    reference leaves behind, so the zero point is taken as a median and the slope is
    measured on what is left.
    """
    zp = reference - m_inst
    offset = np.median(zp)
    resid = zp - offset
    keep = np.abs(resid - np.median(resid)) < 3 * np.std(resid)
    slope, _ = np.polyfit(colour[keep], resid[keep], 1), None
    err = np.sqrt(np.diag(np.polyfit(colour[keep], resid[keep], 1, cov=True)[1]))[0]
    return {"zp": float(offset), "resid": resid, "keep": keep,
            "slope": float(slope[0]), "slope_err": float(err),
            "rms": float(np.std(resid[keep])), "n": int(keep.sum()),
            # Unclipped too, because that is what bayer_zp_test.ecsv records; quoting
            # only the clipped number would make this frame look better than that
            # survey by a definition rather than by a measurement.
            "rms_raw": float(np.std(resid)), "n_raw": int(resid.size)}


# --- figures ---------------------------------------------------------------------------------

#: The two calibration schemes being compared, everywhere they appear. Slots 0 and 1 of
#: the house categorical palette, validated as a pair (worst adjacent CVD dE 24.7,
#: protan). Band identity is carried by the panel title, never by these hues, so the
#: same two colours mean the same two things in every panel.
SCHEME = {
    "jkc": {"color": _style.CATEGORICAL[1], "marker": "o", "ls": "--",
            "label": "classical reference (JKC)"},
    "s50": {"color": _style.CATEGORICAL[0], "marker": "s", "ls": "-",
            "label": "S50 reference"},
}

#: Single-hue light-to-dark ramp for the density panels, stepped off the house blue.
DENSITY = ("#e8f0fa", "#a9c9ee", "#5f9ae0", "#2a78d6", "#17457c")


def _density_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("s50_density", DENSITY)


def figure_passbands(out_dir):
    """Where the S50's bands sit against the classical and Sloan systems."""
    import matplotlib.pyplot as plt
    from astropy.table import Table

    responses = s50.load_passbands()
    standard = Table.read(str(STANDARD))
    wl = s50.XP_GRID

    panels = [("Johnson-Kron-Cousins", ("B_JKC", "V_JKC", "R_JKC", "I_JKC")),
              ("SDSS", ("g_SDSS", "r_SDSS", "i_SDSS"))]
    with plt.rc_context(_style.rc()):
        fig, axes = plt.subplots(2, 1, figsize=(7.4, 5.6), sharex=True)
        for ax, (title, bands) in zip(axes, panels):
            for band in S50_BANDS:
                ax.fill_between(wl, responses[band], color=_style.BAND_COLOR[band],
                                alpha=0.18, lw=0)
                ax.plot(wl, responses[band], color=_style.BAND_COLOR[band],
                        ls=_style.BAND_STYLE[band], lw=1.8,
                        label=f"{band}$_{{S50}}$" if ax is axes[0] else None)
            for band in bands:
                curve = np.asarray(standard[band], float)
                ax.plot(wl, curve, color=_style.REFERENCE, lw=1.0, alpha=0.85)
                peak = wl[int(np.argmax(curve))]
                ax.text(peak, curve.max() + 0.04, band.split("_")[0],
                        ha="center", va="bottom", fontsize=7.5,
                        color=_style.INK_SECONDARY)
            ax.set_ylim(0, 1.22)
            ax.set_ylabel("relative response")
            ax.set_title(f"S50 against {title}", loc="left")
        axes[-1].set_xlim(340, 1020)
        axes[-1].set_xlabel("wavelength [nm]")
        _style.legend(axes[0], ncol=3)
        fig.tight_layout()
        fig.savefig(out_dir / "passbands.png")
        plt.close(fig)


def figure_transformations(out_dir, sample, fits):
    """Each conversion as a density, with the linear and cubic fits over it."""
    import matplotlib.pyplot as plt

    colour = (np.asarray(sample[f"{COLOUR[0]}_S50"])
              - np.asarray(sample[f"{COLOUR[1]}_S50"]))
    grid = np.linspace(*COLOUR_RANGE, 200)
    by_key = {(r["s50_band"], r["target"], r["degree"]): r for r in fits}

    with plt.rc_context(_style.rc()):
        fig, axes = plt.subplots(2, 3, figsize=(9.6, 5.4))
        for row, system in enumerate(("JKC", "SDSS")):
            for col, band in enumerate(("B", "G", "R")):
                target = [t for t in PARTNER[band] if t.endswith(system)][0]
                ax = axes[row, col]
                delta = np.asarray(sample[target]) - np.asarray(sample[f"{band}_S50"])
                good = (np.isfinite(delta) & (colour >= COLOUR_RANGE[0])
                        & (colour <= COLOUR_RANGE[1]))
                lo, hi = np.percentile(delta[good], [0.2, 99.8])
                pad = 0.08 * (hi - lo)
                ax.hexbin(colour[good], delta[good], gridsize=90, bins="log",
                          cmap=_density_cmap(), linewidths=0, mincnt=1,
                          extent=(*COLOUR_RANGE, lo, hi))
                for degree, style in ((1, "--"), (3, "-")):
                    row_fit = by_key[(band, target, degree)]
                    ax.plot(grid, np.polyval(row_fit["coef"], grid), ls=style, lw=1.6,
                            color=_style.MODEL,
                            label=f"degree {degree}" if (row == 0 and col == 0) else None)
                ax.set_xlim(*COLOUR_RANGE)
                ax.set_ylim(lo - pad, hi + pad)
                d1 = 1000 * by_key[(band, target, 1)]["rms"]
                d3 = 1000 * by_key[(band, target, 3)]["rms"]
                # The stats block goes in whichever bottom corner the trend leaves
                # empty, which the sign of the linear term already tells us.
                rising = by_key[(band, target, 1)]["coef"][0] > 0
                _style.annotate(ax, f"linear  {d1:.0f} mmag\ncubic   {d3:.0f} mmag",
                                loc="lower right" if rising else "lower left")
                ax.set_title(f"{target.split('_')[0]} $-$ {band}$_{{S50}}$", loc="left")
                if row == 1:
                    ax.set_xlabel("$B_{S50} - R_{S50}$")
                if col == 0:
                    ax.set_ylabel(f"{system} $-$ S50 [mag]")
        _style.legend(axes[0, 0])
        fig.tight_layout()
        fig.savefig(out_dir / "transformations.png")
        plt.close(fig)


def figure_polynomial_order(out_dir, fits):
    """Does the conversion need more than a linear colour term?"""
    import matplotlib.pyplot as plt

    with plt.rc_context(_style.rc()):
        fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.4), sharey=True)
        for ax, system in zip(axes, ("JKC", "SDSS")):
            for slot, band in enumerate(("B", "G", "R")):
                target = [t for t in PARTNER[band] if t.endswith(system)][0]
                rms = [1000 * r["rms"] for d in DEGREES for r in fits
                       if (r["s50_band"], r["target"], r["degree"]) == (band, target, d)]
                ax.plot(DEGREES, rms, color=_style.BAND_COLOR[band],
                        ls=_style.BAND_STYLE[band], marker=_style.BAND_MARKER[band],
                        ms=6, label=f"{band}$_{{S50}}$ $\\to$ {target.split('_')[0]}")
            ax.set_xticks(DEGREES)
            ax.set_xlabel("polynomial degree in $B_{S50} - R_{S50}$")
            # Log, because the bands span a factor of ten and the question is "by what
            # factor does another term help" -- but with plain tick labels, since
            # matplotlib's default 4x10^1 is unreadable at this size.
            ax.set_yscale("log")
            ax.set_yticks([5, 10, 20, 50])
            ax.set_yticklabels(["5", "10", "20", "50"])
            ax.minorticks_off()
            ax.set_title(system, loc="left")
            # Not the house helper: its 2.5x marker scale exists for dense translucent
            # scatters and is oversized on a line chart whose markers are already 6 pt.
            ax.legend(markerscale=1.0)
        axes[0].set_ylabel("residual rms [mmag]")
        _style.reference_line(axes[1], 10.0, label="10 mmag")
        fig.tight_layout()
        fig.savefig(out_dir / "polynomial-order.png")
        plt.close(fig)


def figure_frame_zeropoint(out_dir, planes):
    """The payoff: residual against colour for one stack, both references."""
    import matplotlib.pyplot as plt

    with plt.rc_context(_style.rc()):
        fig, axes = plt.subplots(1, 3, figsize=(9.8, 3.6), sharex=True, sharey=True)
        for ax, plane in zip(axes, ("B", "G", "R")):
            data = planes[plane]
            _style.reference_line(ax, 0.0)
            stats = []
            for kind, ref in zip(("jkc", "s50"), REFERENCES[plane]):
                fit = zeropoint_fit(data[ref], data["m_inst"], data["bp_rp"])
                style = SCHEME[kind]
                keep = fit["keep"]
                ax.plot(data["bp_rp"][keep], fit["resid"][keep], style["marker"],
                        ms=3.4, alpha=0.5, color=style["color"], mew=0,
                        label=style["label"] if plane == "B" else None)
                grid = np.linspace(data["bp_rp"][keep].min(),
                                   data["bp_rp"][keep].max(), 50)
                ax.plot(grid, fit["slope"] * (grid - np.median(data["bp_rp"][keep])),
                        ls=style["ls"], lw=1.8, color=style["color"])
                stats.append(f"{'JKC' if kind == 'jkc' else 'S50'}  "
                             f"{1000*fit['rms']:3.0f} mmag  "
                             f"{fit['slope']:+.3f} mag/mag")
            ax.set_title(f"{plane} plane", loc="left")
            ax.set_xlabel("$BP - RP$")
            _style.annotate(ax, "\n".join(stats), loc="lower right")
        axes[0].set_ylabel("zero-point residual [mag]")
        axes[0].set_ylim(-0.95, 0.95)
        # One legend for the whole row, above the figure rather than above an axes: the
        # two schemes mean the same thing in every panel, and an axes-level outside
        # legend lands in the strip the panel title occupies.
        handles, labels = axes[0].get_legend_handles_labels()
        fig.tight_layout(rect=(0, 0, 1, 0.91))
        fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
                   fontsize=8, markerscale=2.2, bbox_to_anchor=(0.5, 1.0))
        fig.savefig(out_dir / "frame-zeropoint.png")
        plt.close(fig)


#: A deep 20 Cep stack for the zero-point-linearity figure. 60 minutes rather than the
#: 15 used elsewhere on purpose: the question is how the zero point behaves *across the
#: magnitude range*, and a 47-minute-deeper stack reaches ~2.5 mag fainter, which is
#: where the two references stop agreeing.
DEEP_FRAME = Path("D:/seestar_saturation_data_set/stacks_local/n2_c17_dither/X60min/"
                  "n2_c17_ditherlocal_X060min_20260813-0000_n141_t2820s_20s_IRCUT.fit")

#: References compared in that figure, both for the green plane: the status quo and the
#: S50's own band. Green because that is the science band. The third entry is a plain-text
#: name for the stats block, which is not a mathtext context.
LINEARITY_REFS = (("jkc", "v_jkc_mag", "Gaia synthetic $V$", "Gaia synthetic V"),
                  ("s50", "g_s50_mag", "$G_{S50}$", "G_S50"))


def figure_zp_linearity(out_dir, data, frame_name):
    """Per-star zero point against reference magnitude, for both green references.

    The single most diagnostic view of a calibration: a flat cloud is a zero point that
    means the same thing at every brightness, and every way it can fail shows up as a
    different shape. Saturation bends the bright end up, the noise floor flares the faint
    end out, and a reference on the wrong system tilts the whole thing.
    """
    import matplotlib.pyplot as plt

    finite = np.isfinite(data["m_inst"])
    cal = calibrator_mask(data)

    with plt.rc_context(_style.rc()):
        fig, axes = plt.subplots(2, 1, figsize=(7.6, 6.4), sharex=True,
                                 gridspec_kw={"height_ratios": [2.0, 1.0]})
        stats = []
        for kind, column, label, plain in LINEARITY_REFS:
            style = SCHEME[kind]
            ref = data[column]
            good = finite & np.isfinite(ref)
            zp = ref - data["m_inst"]
            fitted = np.median(zp[cal])
            resid = zp[cal] - fitted
            steady = np.abs(resid) < 3 * np.std(resid)
            rms = float(np.std(resid[steady]))

            axes[0].plot(ref[good], zp[good], style["marker"], ms=3.0, alpha=0.35,
                         color=style["color"], mew=0, label=label)
            _style.reference_line(axes[0], fitted, style=style["ls"])
            anchor = (fitted, ref[good], zp[good])
            axes[1].plot(ref[cal], resid, style["marker"], ms=3.4, alpha=0.5,
                         color=style["color"], mew=0)
            stats.append(f"{plain}:  ZP {fitted:.3f}   {1000*rms:.0f} mmag "
                         f"over {int(steady.sum())} stars")

        # The calibrator window, which is what actually protects the fit from both ends.
        for ax in axes:
            ax.axvspan(*CAL_RANGE, color=_style.GRID, alpha=0.55, lw=0, zorder=0)

        axes[0].set_ylabel("per-star zero point   ref $-$ $m_{inst}$ [mag]")
        axes[0].set_title(f"20 Cep, {frame_name}", loc="left")
        _style.annotate(axes[0], "\n".join(stats), loc="lower right")
        _style.legend(axes[0], ncol=2, loc="upper left")

        # No saturation line: on a stack of 141 subs the peak tops out at 53140 ADU, so
        # no ADU test fires, yet the flux was already clipped in the individual subs.
        # The bright end is labelled for what it is rather than by a threshold that
        # cannot see it.
        fitted, ref, zp = anchor
        # Wide enough to show both failure modes, with the detail left to the lower
        # panel: the flat region is 0.05 mag tall and would be a line at this scale.
        axes[0].set_ylim(fitted - 2.7, fitted + 2.7)
        left = float(np.nanpercentile(ref, 0.2)) - 0.4
        axes[0].set_xlim(left, float(np.nanpercentile(ref, 99.9)) + 0.4)
        # Directly under the bright cloud rather than on an arrow: the points it names
        # are the only thing out there, so a leader line is clutter.
        axes[0].text(left + 0.1, fitted - 1.05, "clipped in the subs",
                     ha="left", va="top", fontsize=7.5, color=_style.INK_SECONDARY)
        axes[1].set_ylim(-0.26, 0.26)
        axes[1].set_ylabel("residual [mag]")
        axes[1].set_xlabel("reference magnitude [mag]")
        _style.reference_line(axes[1], 0.0)
        _style.annotate(axes[1], "calibrators only, each about its own zero point",
                        loc="upper left")
        fig.tight_layout()
        fig.savefig(out_dir / "zp-linearity.png")
        plt.close(fig)


def figure_survey(out_dir):
    """The same comparison over all 57 stacks, so one frame is not the whole claim."""
    import matplotlib.pyplot as plt
    from astropy.table import Table

    path = Path(r"D:/CrowdSky/crowdsky-workspace/Repos/photometry-adventures-2"
                r"/paper/data/bayer_zp_test.ecsv")
    survey = Table.read(str(path))
    rng = np.random.default_rng(3)

    with plt.rc_context(_style.rc()):
        fig, ax = plt.subplots(figsize=(6.6, 3.6))
        for i, plane in enumerate(("B", "G", "R")):
            for j, kind in enumerate(("jkc", "s50")):
                name = {"jkc": "jkc", "s50": "bayer"}[kind]
                rows = survey[(survey["plane"] == plane)
                              & (survey["reference"] == name)]
                values = 1000 * np.asarray(rows["rms"], float)
                y = i + (j - 0.5) * 0.34 + rng.normal(0, 0.035, len(values))
                style = SCHEME[kind]
                ax.plot(values, y, style["marker"], ms=3.6, alpha=0.45,
                        color=style["color"], mew=0,
                        label=style["label"] if i == 0 else None)
                ax.plot(np.median(values), i + (j - 0.5) * 0.34, "|", ms=16, mew=2.2,
                        color=style["color"])
        ax.set_yticks(range(3))
        ax.set_yticklabels(["B plane", "G plane", "R plane"])
        ax.set_xlabel("calibrator rms of the zero point [mmag]")
        ax.set_title("57 fifteen-minute 20 Cep stacks", loc="left")
        _style.legend(ax)
        fig.tight_layout()
        fig.savefig(out_dir / "survey-scatter.png")
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--src", type=Path, required=True, help="the v2 dataset")
    ap.add_argument("--files", type=int, default=24, help="bulk files to sample")
    ap.add_argument("--refresh", action="store_true", help="rebuild the cached sample")
    args = ap.parse_args()

    from astropy.table import Table

    cache = args.work / "s50_system.ecsv"
    if cache.exists() and not args.refresh:
        sample = Table.read(str(cache))
        print(f"sample: {len(sample)} stars (cached)")
    else:
        sample = build_sample(args.work, args.src, args.files)
        sample.write(str(cache), overwrite=True)

    offset, scatter, n = gaia_scale_offset(sample)
    print(f"\nsynthetic V vs Gaia GSPC v_jkc_mag: offset {offset:+.4f} mag "
          f"(AB-Vega), scatter {1000*scatter:.1f} mmag over {n} stars")

    print("\ntransformations, delta = target - S50, in (B_S50 - R_S50):")
    for row in transformations(sample):
        terms = "  ".join(f"{c:+.5f}" for c in row["coef"])
        print(f"  {row['s50_band']}_S50 -> {row['target']:7s} deg {row['degree']}  "
              f"rms {1000*row['rms']:6.1f} mmag  n {row['n']:6d}  [{terms}]")

    print(f"\nzero point of {FRAME.name}:")
    planes = frame_calibrators(FRAME, args.src)
    for plane, data in planes.items():
        for kind, ref in zip(("JKC ", "S50 "), REFERENCES[plane]):
            fit = zeropoint_fit(data[ref], data["m_inst"], data["bp_rp"])
            print(f"  {plane} vs {kind}({ref:11s}) n {fit['n']:4d}  ZP {fit['zp']:7.3f}  "
                  f"rms {1000*fit['rms']:6.1f} mmag (raw {1000*fit['rms_raw']:6.1f})  "
                  f"slope {fit['slope']:+.4f} +/- {fit['slope_err']:.4f} mag/mag")

    print("\nband effective wavelengths and widths [nm]:")
    for row in band_statistics():
        print(f"  {row['band']:8s} lam_eff {row['lam_eff']:6.1f}  FWHM {row['fwhm']:5.1f}")

    FIGURES.mkdir(parents=True, exist_ok=True)
    fits = transformations(sample)
    write_coefficients(FIGURES.parent, fits)
    figure_passbands(FIGURES)
    figure_transformations(FIGURES, sample, fits)
    figure_polynomial_order(FIGURES, fits)
    figure_frame_zeropoint(FIGURES, planes)

    print(f"\nzero-point linearity on {DEEP_FRAME.name}:")
    deep = measure_frame(DEEP_FRAME, args.src)["G"]
    figure_zp_linearity(FIGURES, deep, DEEP_FRAME.name)
    figure_survey(FIGURES)
    print(f"\nfigures -> {FIGURES}")


if __name__ == "__main__":
    main()
