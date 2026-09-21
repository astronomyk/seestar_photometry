#!/usr/bin/env python
"""Fit each Seestar model's total-system passbands from its own frames.

The sensor QE curves (``tools/digitise_sensor_qe.py``) are only half of a passband. The
other half is the UV/IR-cut filter in front of it, and its red edge is not published for
any of these units -- so it is **fitted from the photometry**, exactly as it was for the
S50: calibrate each Bayer plane against a synthetic magnitude computed at a trial edge,
measure the zero-point residual's slope against colour, and take the edge where that
slope crosses zero. A passband that is too red makes the plane look too sensitive to red
stars and the slope goes one way; too blue and it goes the other.

The colour axis is Gaia ``BP-RP``, deliberately *not* a colour built from the bands being
fitted, which would be scored against itself.

Two free parameters per unit, fitted in that order:

``ircut_red_edge``
    Scanned per plane. G and R pin it -- their passbands sit under the edge and move
    most -- and the adopted value is their mean.
``b_tail_scale``
    With the edge fixed, the B plane's residual slope is scanned against a multiplicative
    boost to the blue dye's red tail. The datasheet plots understate that tail (it is the
    thinnest part of a printed curve), and on the S50 the B plane needed +15%.

Two stages, cached, because the first is minutes and the second is milliseconds:

``measure``
    Solve and do forced photometry on every stack of every unit; write per-star
    instrumental magnitudes to ``<work>/photometry_<unit>.ecsv``. Independent of any
    passband, which is the point -- the scan reuses one measurement.
``scan``
    Integrate the field's XP spectra through each trial passband and find the crossings.

Usage::

    uv run python tools/fit_unit_passbands.py --work D:/tmp/xp_s50 --stage all

Maintenance script, not shipped in the wheel.
"""

import argparse
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import build_gaia_s50 as s50  # noqa: E402
from seestar_photometry import _style  # noqa: E402

HERE = Path(__file__).resolve().parent
SENSOR_QE = HERE / "sensor_qe.ecsv"
DATA = Path("D:/seestar_saturation_data_set")


@dataclass
class Unit:
    """One telescope model: its sensor, its data, and what its passbands are called."""

    name: str
    sensor: str
    #: Dataset directories to take 20 Cep stacks from. More than one is more than one
    #: physical unit of the same model, which is the point -- a parameter fitted on a
    #: single serial number is a property of that serial number until shown otherwise.
    datasets: tuple = ()
    subdir: str = "20 Cep"
    #: Explicit directories, for data that does not sit under ``DATA/<dataset>/<subdir>``.
    roots: tuple = ()

    def paths(self):
        """Every stack for this unit, in a stable order."""
        here = [DATA / d / self.subdir for d in self.datasets] + [Path(r) for r in
                                                                  self.roots]
        return [sorted(d.glob("*.fit")) for d in here]


UNITS = (
    # The control. Its passbands are already published and were fitted by a different
    # script in a different repository, so running it through *this* code path is the
    # check that the method reproduces a known answer before it is believed on a sensor
    # where there is nothing to compare against. The same 15-minute local stacks the
    # published fit used.
    Unit("s50", "imx462", roots=(
        DATA / "stacks_local" / "n2_c17_dither" / "X15min",
        DATA / "stacks_local" / "n2_b61_altaz" / "X15min",
    )),
    Unit("s30pro", "imx585", ("s30pro_77febe82", "s30pro_ee65a72a")),
    Unit("s50pro", "os08b10", ("s50pro_39ac5218", "s50pro_63ea2f59", "s50pro_ded957a0")),
)

#: The S50's fitted model, carried here so the new units are fitted the same way and can
#: be compared against it. See ``photometry-adventures-2/paper/scripts/``.
IRCUT_BLUE_EDGE = 390.0
IRCUT_EDGE_WIDTH = 15.0
IRCUT_PEAK = 0.97
#: Where each channel's "red tail" correction ramps in: flat below, full scale above.
#: Placed just past each dye's own peak, so the correction touches only the weak shoulder
#: and leaves the band's core alone.
TAIL_RAMP = {"G": (580.0, 620.0), "B": (520.0, 540.0)}

#: Search grids. The red-edge range brackets every plausible IR-cut; the tail range runs
#: from "the datasheet is right" upward, since a printed curve understates a thin tail
#: rather than overstating it.
EDGE_GRID = np.arange(620.0, 800.0, 1.0)
TAIL_GRID = np.arange(1.00, 6.01, 0.02)

#: Calibrator cuts, matching ``build_s50_system_figures.py`` so the numbers are readable
#: against the S50's.
CAL_RANGE = (9.5, 13.0)
SNR_MIN = 30.0
SAT_FRAC = 0.90
CLIP_ADU = 65535.0
ENCLOSED = 0.95
BLEND_APER = 1.5
BLEND_DMAG = 2.5

#: The fitted result, pasted back in once the scan has run, so the passband file can be
#: rebuilt without repeating two hours of photometry. Every number here came out of
#: ``--stage scan`` on the datasets in :data:`UNITS`.
#:
#: The ``s50`` row is the **control** and is not what the package ships: the published
#: S50 passbands (``tools/bayer_passbands.ecsv``, edge 693.0, B tail 1.15, no G tail)
#: are what the 34.2M synthetic magnitudes in the catalogue were built from, and
#: rebuilding them from a refit would silently invalidate those. It is here because
#: reproducing it -- 692.8 nm against 693.0, and a G tail of 1.016 where the published
#: model has none -- is what says the method still works.
FITTED = {
    "s50": {"red_edge": 692.8, "tails": {"G": 1.016, "B": 1.204}},
    "s30pro": {"red_edge": 678.4, "tails": {"G": 1.723, "B": 2.452}},
    "s50pro": {"red_edge": 684.9, "tails": {"G": 1.644, "B": 2.043}},
}

#: Field centre and radius for the Gaia cone, in degrees. Covers the S30 Pro's 2.2 x 4.0
#: degree field, which is the largest of the three.
FIELD_CENTRE = (331.53, 62.92)
FIELD_RADIUS = 2.4


# --- the passband model ---------------------------------------------------------------

def ircut(wl, red_edge):
    """Logistic-edged flat top, scaled so the 10-90% width is :data:`IRCUT_EDGE_WIDTH`."""
    k = 2 * np.log(9.0) / IRCUT_EDGE_WIDTH
    rise = 1.0 / (1.0 + np.exp(-k * (wl - IRCUT_BLUE_EDGE)))
    fall = 1.0 / (1.0 + np.exp(k * (wl - red_edge)))
    return IRCUT_PEAK * rise * fall


def passbands(qe, wl, red_edge, tails=None):
    """``{band: response}`` for one unit at a trial edge and set of tail scales.

    ``qe`` is ``{band: sensor QE}`` on ``wl``; ``tails`` is ``{band: scale}`` for the
    red-shoulder corrections, which default to 1 (the datasheet as drawn).
    """
    cut = ircut(wl, red_edge)
    tails = tails or {}
    out = {}
    for band in ("R", "G", "B"):
        response = np.asarray(qe[band], float).copy()
        scale = float(tails.get(band, 1.0))
        if scale != 1.0 and band in TAIL_RAMP:
            lo, hi = TAIL_RAMP[band]
            ramp = np.clip((wl - lo) / (hi - lo), 0.0, 1.0)
            response = response * (1.0 + (scale - 1.0) * ramp)
        out[band] = response * cut
    peak = max(v.max() for v in out.values())
    return {band: v / peak for band, v in out.items()}


def sensor_qe(sensor, wl):
    """One sensor's digitised QE, resampled onto ``wl``."""
    from astropy.table import Table

    table = Table.read(str(SENSOR_QE))
    source = np.asarray(table["wavelength_nm"], float)
    return {band: np.interp(wl, source, np.asarray(table[f"{band}_{sensor}"], float))
            for band in ("R", "G", "B")}


# --- stage 1: photometry ----------------------------------------------------------------

def field_catalogue(gaiadb_dir, epoch=None):
    """Gaia sources over the field, from the **offline** catalogue.

    Not a TAP cone. A 4.6 degree box at G < 15 is nine tiles and the endpoint was still
    chewing on it after ten minutes; the offline copy answers the same question from
    disk in seconds, which is the entire reason it exists. It also carries proper
    motions, so the positions can be moved to the observing epoch -- Gaia is J2016.0 and
    these frames are 2026.7, which is an arcsecond for a 100 mas/yr star.

    Two columns it cannot supply, because they live only in ``gaia_source``:

    ``bp_rp``
        Synthesised instead, from the same XP spectra as everything else, through Gaia's
        own BP and RP passbands. See :func:`add_synthetic_colour`, which checks it
        against real Gaia colours before it is used.
    ``phot_variable_flag``
        Not recoverable offline, so the variable cut is simply not applied. It removes a
        handful of stars from a few hundred; the 3-sigma clip in the slope fit is what
        actually protects against the one that matters.
    """
    from seestar_photometry import gaiadb

    return gaiadb.cone(FIELD_CENTRE, FIELD_RADIUS, directory=gaiadb_dir,
                       gmag_limit=15.0, epoch=epoch)


def add_synthetic_colour(catalogue, ids, flux, flux_err, quiet=False):
    """Attach a ``bp_rp`` synthesised from the XP spectra, and check it against Gaia's.

    The check is not optional decoration: a colour axis that is quietly wrong would tilt
    every fitted edge in the same direction and still look perfectly smooth. The 20 Cep
    TAP result in ``photometry-adventures-2`` carries real Gaia ``bp_rp`` for this field,
    so the synthetic colour is compared against it star by star and the agreement is
    reported.
    """
    from astropy.table import Table

    standard = Table.read(str(HERE / "standard_passbands.ecsv"))
    colour = {}
    for band in ("BP_GAIA", "RP_GAIA"):
        response = np.asarray(standard[band], float)
        colour[band], _ = s50.synth_mags(flux, flux_err, s50.XP_GRID, response)
    synthetic = colour["BP_GAIA"] - colour["RP_GAIA"]

    star = np.asarray(catalogue["source_id"], np.int64)
    where = np.clip(np.searchsorted(ids, star), 0, max(len(ids) - 1, 0))
    hit = len(ids) and (ids[where] == star)
    out = np.full(len(catalogue), np.nan)
    out[hit] = synthetic[where[hit]]
    catalogue["bp_rp"] = out

    truth = Path("D:/CrowdSky/crowdsky-workspace/Repos/photometry-adventures-2"
                 "/paper/data/gaia_20Cep.ecsv")
    if truth.exists() and not quiet:
        real = Table.read(str(truth))["source_id", "bp_rp"]
        lut = {int(s): float(c) for s, c in zip(real["source_id"], real["bp_rp"])}
        common = np.array([lut.get(int(s), np.nan) for s in star])
        good = np.isfinite(common) & np.isfinite(out)
        if good.sum() > 100:
            delta = out[good] - common[good]
            print(f"  synthetic BP-RP vs Gaia's own, {int(good.sum())} stars: "
                  f"median {np.median(delta):+.4f}, scatter "
                  f"{np.std(delta - np.median(delta)):.4f} mag", flush=True)
    return catalogue


def _blend_flag(inframe, catalogue, aperture_px, pixel_scale):
    """True where a catalogue neighbour is close and bright enough to pollute the sum."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    radius = BLEND_APER * aperture_px * pixel_scale
    here = SkyCoord(inframe["ra"], inframe["dec"], unit="deg")
    everything = SkyCoord(catalogue["ra"], catalogue["dec"], unit="deg")
    g_all = np.asarray(catalogue["phot_g_mean_mag"], float)
    g_here = np.asarray(inframe["phot_g_mean_mag"], float)
    i_all, i_here, sep, _ = here.search_around_sky(everything, radius * u.arcsec)
    flag = np.zeros(len(inframe), bool)
    real = sep.arcsec > 1e-6
    contaminating = g_all[i_all[real]] < g_here[i_here[real]] + BLEND_DMAG
    flag[i_here[real][contaminating]] = True
    return flag


def measure_unit(unit, catalogue, limit=0, quiet=False):
    """Per-star instrumental magnitudes for every stack of one unit."""
    from astropy.table import Table, vstack
    from seestar_photometry import astrometry, catalogs, frames, photometry

    rows = []
    for paths in unit.paths():
        if limit:
            paths = paths[:limit]
        for path in paths:
            try:
                frame = frames.load_frame(str(path))
                wcs = astrometry.solve(frame, solver="astap")
            except Exception as exc:                       # noqa: BLE001
                print(f"  {path.name}: load/solve failed ({exc})", flush=True)
                continue
            if wcs is None:
                print(f"  {path.name}: no WCS", flush=True)
                continue
            inframe = catalogs.sources_in_frame(catalogue, wcs, frame)
            # Only stars that could ever be calibrators are measured. The cut window is
            # V in [9.5, 13); everything fainter is 90% of the positions on a 4K frame
            # and none of the answer, and forced photometry is linear in positions.
            # Blends are still flagged against the *whole* catalogue below -- a faint
            # neighbour still pollutes an aperture even though it is never a calibrator.
            bright = np.ma.filled(np.asarray(inframe["v_jkc_mag"]), np.inf) \
                < CAL_RANGE[1] + 0.5
            inframe = inframe[bright]
            if len(inframe) < 50:
                print(f"  {path.name}: only {len(inframe)} candidate calibrators",
                      flush=True)
                continue
            phot = photometry.forced_photometry(
                frame, inframe["ra"], inframe["dec"], wcs,
                source_id=inframe["source_id"], enclosed=ENCLOSED)
            blend = _blend_flag(inframe, catalogue, phot.meta["aperture"][1],
                                astrometry.pixel_scale(frame))
            lut = {int(s): i for i, s in enumerate(inframe["source_id"])}

            for plane in ("R", "G", "B"):
                sub = phot[(phot["band"] == plane) & phot["on_chip"]]
                index = np.array([lut[int(s)] for s in sub["source_id"]])
                flux = np.asarray(sub["flux"], float)
                with np.errstate(invalid="ignore", divide="ignore"):
                    m_inst = -2.5 * np.log10(np.where(flux > 0, flux, np.nan))

                def column(name, index=index):
                    return np.ma.filled(
                        np.asarray(inframe[name])[index], np.nan).astype(float)

                rows.append(Table({
                    "source_id": np.asarray(sub["source_id"], np.int64),
                    "frame": np.full(len(sub), path.name),
                    "unit": np.full(len(sub), unit.name),
                    "plane": np.full(len(sub), plane),
                    "m_inst": m_inst,
                    "bp_rp": column("bp_rp"),
                    "v_jkc_mag": column("v_jkc_mag"),
                    "snr": np.asarray(sub["snr"], float),
                    "peak": np.asarray(sub["max_pix_value"], float),
                    "blend": blend[index],
                }))
            if not quiet:
                print(f"  {path.name}: {len(inframe)} sources, "
                      f"aperture {phot.meta['aperture'][1]:.2f} px", flush=True)
    if not rows:
        raise RuntimeError(f"{unit.name}: no frames measured")
    return vstack(rows)


def calibrator_mask(table):
    """Bright enough, unsaturated, unblended, steady, and with a colour to fit against."""
    return (np.isfinite(np.asarray(table["m_inst"], float))
            & np.isfinite(np.asarray(table["bp_rp"], float))
            & (np.asarray(table["snr"], float) > SNR_MIN)
            & (np.asarray(table["peak"], float) < SAT_FRAC * CLIP_ADU)
            & (np.asarray(table["v_jkc_mag"], float) >= CAL_RANGE[0])
            & (np.asarray(table["v_jkc_mag"], float) < CAL_RANGE[1])
            & ~np.asarray(table["blend"], bool))


# --- stage 2: the scan --------------------------------------------------------------------

def load_spectra(work, source_ids, quiet=False):
    """``(ids, flux, flux_err)`` for the field, from the local bulk XP files."""
    wanted = np.unique(np.asarray(source_ids, np.int64))
    entries = s50.enumerate_files(work)
    hp8 = np.unique(wanted >> s50.HPX8_SHIFT)
    needed = [row for row in entries if np.any((hp8 >= row[1]) & (hp8 <= row[2]))]
    ids, flux, error = [], [], []
    for name, _, _, size in needed:
        path = s50.local_path(name, Path(work) / "files")
        if not path.exists():
            path, _ = s50.fetch_one(name, size, Path(work) / "files", quiet)
        got, f, e = s50.parse_file(path, wanted)
        if len(got):
            ids.append(got)
            flux.append(f)
            error.append(e)
    if not ids:
        raise RuntimeError("no XP spectra for this field")
    ids = np.concatenate(ids)
    order = np.argsort(ids)
    if not quiet:
        print(f"  {len(ids)} XP spectra from {len(needed)} bulk files", flush=True)
    return ids[order], np.vstack(flux)[order], np.vstack(error)[order]


def slope_at(table, keep, reference):
    """Zero-point residual slope against BP-RP, after a per-frame median zero point.

    Per *frame*, because the zero point moves with transparency and airmass across a
    night and a single median over all frames would fold that drift into the residual
    and flatten the very slope being measured.
    """
    colour = np.asarray(table["bp_rp"], float)[keep]
    m_inst = np.asarray(table["m_inst"], float)[keep]
    frame = np.asarray(table["frame"])[keep]
    resid = reference - m_inst
    for name in np.unique(frame):
        here = frame == name
        resid[here] -= np.median(resid[here])
    good = np.isfinite(resid) & np.isfinite(colour)
    if good.sum() < 50:
        return np.nan, np.nan, 0
    clip = np.abs(resid - np.median(resid[good])) < 3 * np.std(resid[good])
    good &= clip
    coef, cov = np.polyfit(colour[good], resid[good], 1, cov=True)
    return float(coef[0]), float(np.sqrt(cov[0, 0])), int(good.sum())


def crossing(x, y):
    """Where a monotonic-ish scan crosses zero, by linear interpolation."""
    y = np.asarray(y, float)
    sign = np.sign(y)
    change = np.flatnonzero(sign[:-1] * sign[1:] < 0)
    if not len(change):
        return np.nan
    i = change[0]
    return float(x[i] - y[i] * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))


def scan_unit(unit, table, spectra, quiet=False):
    """Scan the red edge per plane, then the B tail. Returns the fitted parameters."""
    ids, flux, _ = spectra
    qe = sensor_qe(unit.sensor, s50.XP_GRID)

    keeps, refs_index = {}, {}
    star = np.asarray(table["source_id"], np.int64)
    where_all = np.clip(np.searchsorted(ids, star), 0, len(ids) - 1)
    for plane in ("R", "G", "B"):
        keep = (np.asarray(table["plane"]) == plane) & calibrator_mask(table)
        keep &= ids[where_all] == star
        keeps[plane] = keep
        refs_index[plane] = where_all[keep]

    # Only the stars that survive a calibrator cut are ever integrated, and each one
    # once however many rows it has. The scan is ~500 trial passbands; over the whole
    # cone that is 28k spectra a trial, over the calibrators it is a couple of thousand.
    used = np.unique(np.concatenate([refs_index[p] for p in "RGB"]))
    compact = flux[used].astype(np.float64)
    remap = np.full(len(ids), -1, np.int64)
    remap[used] = np.arange(len(used))
    refs_index = {p: remap[refs_index[p]] for p in "RGB"}

    # The magnitude only; synth_mags also propagates the flux errors, which the scan
    # never looks at and which cost as much again.
    f_ab = s50.F_AB * s50.C_NM / s50.XP_GRID**2
    reference_cache = {}

    def reference_for(plane, red_edge, tails=None):
        key = (plane, round(red_edge, 4),
               tuple(sorted((k, round(v, 4)) for k, v in (tails or {}).items())))
        if key not in reference_cache:
            response = passbands(qe, s50.XP_GRID, red_edge, tails)[plane]
            weight = response * s50.XP_GRID
            signal = np.trapezoid(compact * weight, s50.XP_GRID, axis=1)
            ref = np.trapezoid(f_ab * weight, s50.XP_GRID)
            with np.errstate(invalid="ignore", divide="ignore"):
                reference_cache[key] = -2.5 * np.log10(signal / ref)
        return reference_cache[key][refs_index[plane]]

    def crossing_over(grid, plane, make):
        slopes = [slope_at(table, keeps[plane], reference_for(*make(value)))[0]
                  for value in grid]
        return crossing(grid, slopes), np.asarray(slopes)

    # The edge comes from R alone. All three planes see the IR cut, but only R's *peak*
    # sits under it -- 70-85% QE right where the filter turns over -- so its crossing is
    # set by the edge and almost nothing else. G and B are in a deep trough there (14-33%
    # and 2-5%), which is the thinnest, steepest part of a printed curve and the easiest
    # few percent to digitise wrong; on both new sensors that is exactly where the model
    # broke. So the edge is measured where the data is strong and the residual is
    # attributed to the shoulders where it is weak, rather than averaged across both and
    # left as an equal-and-opposite colour term in G and R.
    edge, edge_scan = crossing_over(EDGE_GRID, "R", lambda e: ("R", e))
    tails = {}
    for plane in ("G", "B"):
        tails[plane], _ = crossing_over(
            TAIL_GRID, plane,
            lambda t, plane=plane: (plane, edge, {**tails, plane: t}))

    solo = {}
    for plane in ("R", "G", "B"):
        solo[plane], _ = crossing_over(EDGE_GRID, plane, lambda e, p=plane: (p, e))

    residual = {}
    for plane in ("R", "G", "B"):
        slope, err, n = slope_at(table, keeps[plane], reference_for(plane, edge, tails))
        residual[plane] = (slope, err, n)

    if not quiet:
        print(f"    edge from R: {edge:.1f} nm    "
              f"(solo crossings G {solo['G']:.1f}, R {solo['R']:.1f}, "
              f"B {solo['B']:.1f} nm)", flush=True)
        print(f"    tail scales: G x{tails['G']:.3f}, B x{tails['B']:.3f}", flush=True)
        for plane in ("R", "G", "B"):
            slope, err, n = residual[plane]
            print(f"      {plane} residual colour term {slope:+.4f} +/- {err:.4f} "
                  f"mag/mag  (n={n})", flush=True)

    return {"unit": unit.name, "sensor": unit.sensor, "red_edge": edge,
            "tails": tails, "solo_edges": solo, "residual": residual,
            "n": {p: int(keeps[p].sum()) for p in "RGB"}}


def write_passbands(out=HERE / "unit_passbands.ecsv"):
    """Write the fitted total-system passbands for the units this package fits itself.

    The S50 is deliberately absent: its passbands are already published in
    ``bayer_passbands.ecsv`` and the deployed catalogue's magnitudes were integrated
    through that exact file. Two files claiming to be the S50's passbands, differing in
    the third decimal, is how a catalogue quietly stops matching its own documentation.
    """
    from astropy.table import Table

    wl = s50.XP_GRID
    table = Table({"wavelength_nm": wl})
    lines = []
    for unit in UNITS:
        if unit.name == "s50":
            continue
        fit = FITTED[unit.name]
        curves = passbands(sensor_qe(unit.sensor, wl), wl,
                           fit["red_edge"], fit["tails"])
        for band, response in curves.items():
            table[f"{band}_{unit.name}"] = response
        lines.append(
            f"{unit.name} ({unit.sensor}): IR-cut red edge {fit['red_edge']:.1f} nm, "
            f"G tail x{fit['tails']['G']:.3f}, B tail x{fit['tails']['B']:.3f}")

    table.meta["comments"] = [
        "Seestar total-system passbands (relative response), per unit model.",
        "Sensor QE from tools/sensor_qe.ecsv times a UV/IR-cut model whose red edge and",
        "per-channel red-shoulder scales are FITTED from each unit's own 20 Cep frames",
        "by zeroing the zero-point colour term against Gaia BP-RP.",
        f"Blue edge {IRCUT_BLUE_EDGE:.0f} nm, 10-90% width {IRCUT_EDGE_WIDTH:.0f} nm, "
        f"peak {IRCUT_PEAK}.",
        *lines,
        "The S50 is NOT here: see tools/bayer_passbands.ecsv, which the deployed",
        "catalogue's r/g/b_s50_mag columns were integrated through.",
        f"Built by tools/{Path(__file__).name}.",
    ]
    table.write(str(out), overwrite=True)
    print(f"wrote {out}")
    for line in lines:
        print(f"  {line}")
    return table


FIGURES = Path(__file__).resolve().parents[1] / "docs" / "background-info" / "figures"


def figure_units(out_dir=FIGURES):
    """The three models' passbands, and the sensor curves they came from."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.table import Table

    wl = s50.XP_GRID
    published = Table.read(str(HERE / "bayer_passbands.ecsv"))
    # nan outside the published table's own 390-1000 nm range rather than the default
    # edge-hold, which drew a flat 0.25 blue pedestal down to 350 nm that exists in no
    # file and in no filter.
    rows = [("S50 (IMX462)", {b: np.interp(wl, published["wavelength_nm"],
                                           published[f"{b}_S50"],
                                           left=np.nan, right=np.nan) for b in "RGB"},
             FITTED["s50"]),
            ("S30 Pro (IMX585)", passbands(sensor_qe("imx585", wl), wl,
                                           **FITTED["s30pro"]), FITTED["s30pro"]),
            ("S50 Pro (OS08B10)", passbands(sensor_qe("os08b10", wl), wl,
                                            **FITTED["s50pro"]), FITTED["s50pro"])]

    with plt.rc_context(_style.rc()):
        fig, axes = plt.subplots(3, 1, figsize=(7.6, 7.2), sharex=True)
        for ax, (title, curves, fit) in zip(axes, rows):
            for band in ("R", "G", "B"):
                ax.fill_between(wl, curves[band], color=_style.BAND_COLOR[band],
                                alpha=0.16, lw=0)
                ax.plot(wl, curves[band], color=_style.BAND_COLOR[band],
                        ls=_style.BAND_STYLE[band], lw=1.7,
                        label=f"{band}" if ax is axes[0] else None)
            _style.reference_line(ax, fit["red_edge"], axis="x",
                                  label=f"{fit['red_edge']:.0f} nm")
            ax.set_ylim(0, 1.15)
            ax.set_ylabel("response")
            ax.set_title(title, loc="left")
        _style.legend(axes[0], ncol=3)
        axes[-1].set_xlim(350, 900)
        axes[-1].set_xlabel("wavelength [nm]")
        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "unit-passbands.png")
        plt.close(fig)
    print(f"wrote {out_dir / 'unit-passbands.png'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", type=Path, required=True,
                    help="scratch holding the bulk XP files, as for build_gaia_s50.py")
    ap.add_argument("--gaiadb", type=Path, default=Path("D:/tmp/gaia_s50"),
                    help="the offline Gaia dataset to read the field from")
    ap.add_argument("--epoch", type=float, default=2026.7,
                    help="observing epoch to propagate catalogue positions to")
    ap.add_argument("--stage", choices=("measure", "scan", "write", "all"),
                    default="all")
    ap.add_argument("--limit", type=int, default=0, help="max stacks per dataset")
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()

    from astropy.table import Table

    if args.stage == "write":
        write_passbands()
        figure_units()
        return

    args.work.mkdir(parents=True, exist_ok=True)
    catalogue = field_catalogue(args.gaiadb, epoch=args.epoch)
    print(f"catalogue: {len(catalogue)} Gaia sources over {FIELD_RADIUS} deg "
          f"at {FIELD_CENTRE}, epoch {args.epoch}", flush=True)
    spectra = load_spectra(args.work, np.asarray(catalogue["source_id"]))
    catalogue = add_synthetic_colour(catalogue, *spectra)

    for unit in UNITS:
        cache = args.work / f"photometry_{unit.name}.ecsv"
        if args.stage in ("measure", "all"):
            if cache.exists() and not args.refresh:
                print(f"{unit.name}: cached photometry", flush=True)
            else:
                print(f"{unit.name}: measuring", flush=True)
                t0 = time.time()
                measure_unit(unit, catalogue, limit=args.limit).write(
                    str(cache), overwrite=True)
                print(f"  done in {(time.time()-t0)/60:.1f} min", flush=True)

    if args.stage in ("scan", "all"):
        for unit in UNITS:
            table = Table.read(str(args.work / f"photometry_{unit.name}.ecsv"))
            print(f"\n{unit.name} ({unit.sensor}): {len(table)} measurements, "
                  f"{len(np.unique(table['frame']))} frames", flush=True)
            scan_unit(unit, table, spectra)


if __name__ == "__main__":
    main()
