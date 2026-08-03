"""Photometric calibration onto the Gaia synthetic V system.

Once a frame's sources carry their cross-matched catalogue columns (see
``Extraction.match_gaia``), calibrating the green band onto Gaia DR3 synthetic
Johnson V is a straight-line fit:

    V = m_inst + ZP + k * (B - R - colour0)

where ``m_inst = -2.5 * log10(flux)`` is the fixed-aperture instrumental magnitude,
``B - R`` the synthetic Johnson colour, ``ZP`` the zero point and ``k`` the colour
term. The aperture correction is a constant (each band is sized to its own 90%
enclosure) and is absorbed into ``ZP``.

``ZP`` is quoted at the field's median colour ``colour0``, which decorrelates it from
the colour term -- otherwise the two trade off against each other and the per-frame
zero points are not comparable. The fit is restricted to a magnitude window
(default 10-14: away from saturation at the bright end and the detection floor at the
faint end), and Gaia-flagged variables and sigma-clipped outliers are excluded.

The green Bayer channel through the Seestar's IRCUT filter is a good but not exact
match to Johnson V; the residual mismatch is what the colour term absorbs.
"""

from dataclasses import dataclass

import numpy as np

#: Default magnitude window for the zero-point fit.
FIT_MAG_RANGE = (10.0, 14.0)

#: Fraction of the 16-bit ceiling above which a stacked pixel counts as saturated.
#: Seestar stacks are (scaled) averages, so the hard 65535 clip of a single sub is
#: softened; 0.95 catches the flat-topped bright stars in practice.
SAT_LEVEL = 0.95 * 65535.0

#: Approximate V-band sky-brightness (mag/arcsec^2) lower bounds per Bortle class.
_BORTLE_SQM = [(21.75, 1), (21.5, 2), (21.25, 3), (20.5, 4),
               (19.5, 5), (18.75, 6), (18.0, 7), (17.0, 8)]


@dataclass
class Calibration:
    """Result of :func:`fit_zeropoint`.

    ``V = m_inst + zeropoint + colour_term * (colour - colour0)``.

    ``rms`` is the master photometric-quality number for a frame. ``chi2_red`` is the
    scatter measured against the photon-noise prediction: ~1 means the residuals are
    pure measurement noise, while >>1 flags systematics (cloud, trailing, a bad
    WCS). In practice even excellent Seestar frames sit at 100-200, which is the
    signature of a ~0.03 mag systematic floor rather than a broken fit.
    """

    band: str
    zeropoint: float
    zeropoint_err: float
    colour_term: float
    colour_term_err: float
    colour0: float
    n_stars: int
    rms: float
    chi2_red: float
    mag_range: tuple = FIT_MAG_RANGE
    colour_label: str = ""

    #: Sigma-clipped residuals and the arrays they were fit from, for the
    #: diagnostic figures. Not part of the numeric result.
    fit: dict = None


def _pick_v_column(t):
    """Name of the reference-V column: prefer synthetic JKC V, else plain V."""
    for name in ("v_jkc_mag", "v_mag"):
        if name in t.colnames:
            return name
    raise KeyError("no reference V column (need 'v_jkc_mag' or 'v_mag')")


def _pick_colour(t):
    """``(colour, mask, label)`` from the best available colour source.

    Cascade: synthetic Johnson ``b_jkc_mag - r_jkc_mag`` (best -- it is the same
    photometric system as the reference V) -> Gaia native ``bp_rp`` -> none. The
    Gaia-V vs Bayer-green colour term is small, so a coarser colour (or none) barely
    moves the zero point; this exists so the fit keeps working across catalogue
    schemas rather than to squeeze out accuracy.
    """
    cols = t.colnames
    if "b_jkc_mag" in cols and "r_jkc_mag" in cols:
        colour = np.asarray(t["b_jkc_mag"], float) - np.asarray(t["r_jkc_mag"], float)
        mask = np.ma.getmaskarray(t["b_jkc_mag"]) | np.ma.getmaskarray(t["r_jkc_mag"])
        return colour, mask, "B-R (JKC)"
    if "bp_rp" in cols:
        return np.asarray(t["bp_rp"], float), np.ma.getmaskarray(t["bp_rp"]), "BP-RP"
    return None, None, "none"


def fit_zeropoint(
    sources, band="G", mag_range=FIT_MAG_RANGE, clip=3.0, iters=2,
    exclude_variable=True,
):
    """Fit the zero point (and, where available, a colour term) against reference V.

    The reference V is Gaia synthetic JKC V (``v_jkc_mag``) when present, else a
    plain ``v_mag``. The colour cascades as in :func:`_pick_colour`; with no colour
    at all the fit reduces to a sigma-clipped zero-point mean with
    ``colour_term = 0``.

    Parameters
    ----------
    sources : astropy.table.Table
        Table augmented by ``Extraction.match_gaia`` -- needs ``flux``, ``snr``, a V
        column, and optionally colour columns and ``phot_variable_flag``.
    band : str
        Band to calibrate. ``"G"`` is the science band.
    mag_range : tuple of float or None
        Reference-V window the fit is restricted to. Zero point and colour term are
        fit over the same range so they stay mutually consistent. ``None`` uses all
        matched stars.
    clip, iters : float, int
        Sigma-clip threshold and number of refit iterations.
    exclude_variable : bool
        Drop catalogue-flagged variables from the fit.

    Returns
    -------
    Calibration
    """
    vcol = _pick_v_column(sources)
    t = sources[np.asarray(sources["band"]) == band]
    colour, colour_mask, colour_label = _pick_colour(t)

    finite = ~np.ma.getmaskarray(t[vcol])
    if colour_mask is not None:
        finite &= ~colour_mask
    if exclude_variable and "phot_variable_flag" in t.colnames:
        flag = np.ma.getdata(t["phot_variable_flag"]).astype(str)
        finite &= flag != "VARIABLE"
    if mag_range is not None:
        v_all = np.asarray(t[vcol], dtype=float)
        finite &= (v_all >= mag_range[0]) & (v_all <= mag_range[1])
    # A non-positive forced flux has no instrumental magnitude; excluding it here
    # keeps the fit clean without polluting it with a nan-handling branch.
    finite &= np.asarray(t["flux"], dtype=float) > 0
    t = t[finite]
    colour = colour[finite] if colour is not None else None

    v = np.asarray(t[vcol], dtype=float)
    m_inst = -2.5 * np.log10(np.asarray(t["flux"], dtype=float))
    sigma = 1.0857 / np.asarray(t["snr"], dtype=float)  # photon error on the residual
    residual = v - m_inst  # = ZP + k*colour

    if colour is not None:
        colour0 = float(np.median(colour))
        x = colour - colour0
        keep = np.ones(len(x), dtype=bool)
        for _ in range(iters + 1):
            (k, zp), cov = np.polyfit(x[keep], residual[keep], 1, cov=True)
            model = zp + k * x
            scatter = np.std((residual - model)[keep])
            keep = np.abs(residual - model) < clip * scatter
        k_err = float(np.sqrt(cov[0, 0]))
        zp_err = float(np.sqrt(cov[1, 1]))
    else:
        colour0, k, k_err = float("nan"), 0.0, 0.0
        keep = np.ones(len(residual), dtype=bool)
        for _ in range(iters + 1):
            zp = float(np.mean(residual[keep]))
            scatter = np.std(residual[keep] - zp)
            keep = np.abs(residual - zp) < clip * scatter
        model = np.full(len(residual), zp)
        zp_err = float(np.std(residual[keep]) / np.sqrt(max(int(keep.sum()), 1)))

    npar = 2 if colour is not None else 1
    dof = max(int(keep.sum()) - npar, 1)
    chi2_red = float(np.sum(((residual - model)[keep] / sigma[keep]) ** 2) / dof)

    return Calibration(
        band=band,
        zeropoint=float(zp),
        zeropoint_err=zp_err,
        colour_term=float(k),
        colour_term_err=k_err,
        colour0=colour0,
        n_stars=int(keep.sum()),
        rms=float(np.std((residual - model)[keep])),
        chi2_red=chi2_red,
        mag_range=tuple(mag_range) if mag_range is not None else None,
        colour_label=colour_label,
        fit={
            "v": v, "colour": colour, "m_inst": m_inst, "residual": residual,
            "model": model, "keep": keep, "snr": np.asarray(t["snr"], dtype=float),
            "x": np.asarray(t["x"], dtype=float), "y": np.asarray(t["y"], dtype=float),
        },
    )


def apply_calibration(m_inst, cal, colour=None):
    """Put an instrumental magnitude on the reference V scale.

    ``V = m_inst + ZP + k * (colour - colour0)``. With no colour for the source the
    colour term is dropped, which is correct at the field's median colour and
    degrades linearly away from it -- fine for a red-ish target, worth a thought for
    a very blue one.
    """
    m_inst = np.asarray(m_inst, dtype=float)
    if colour is None or not np.isfinite(cal.colour0):
        return m_inst + cal.zeropoint
    return m_inst + cal.zeropoint + cal.colour_term * (
        np.asarray(colour, dtype=float) - cal.colour0
    )


def limiting_mag(zeropoint, sigma_aper, nsigma=5.0):
    """n-sigma limiting magnitude from the zero point and aperture noise.

    For a background-limited source the aperture noise ``sigma_aper`` (the constant
    aperture flux error, = sqrt(aperture area) x sky RMS) sets the faint limit: the
    source reaches ``nsigma`` when its flux is ``nsigma * sigma_aper``. Quoted at the
    zero point's fiducial colour, where the colour term vanishes.

    Needs no detections near the limit -- it is a noise-model extrapolation of the
    bright-star zero point.

    **Optimistic in absolute terms by ~0.7 mag.** ``sigma_aper`` assumes independent
    pixels, but demosaic interpolation and sub-pixel resampling correlate neighbours:
    the measured aperture noise is ~2.0x the independent-pixel prediction on S50
    on-board stacks. Relative use is fine (the offset is near-constant within one
    processing chain, so exposure-time fits and per-unit rankings stand); absolute
    depth claims are not, and the value must never be compared across processing
    chains. See ``docs/frame-table.md``.
    """
    return float(zeropoint - 2.5 * np.log10(nsigma * sigma_aper))


def saturation_mag(sources, sat_level=SAT_LEVEL):
    """Bright-end (saturation) limit in reference V from a band's matched sources.

    As stars brighten their peak pixel hits the sensor ceiling and the measured flux
    stops tracking V -- the "turn-off" in the V vs instrumental-magnitude relation.
    Located empirically as the *faintest* (largest V) matched source whose peak pixel
    is at or above ``sat_level``: everything brighter is unreliable. Returns nan if
    no matched star saturates, i.e. the limit is brighter than anything measured here.

    ``sources`` must be single-band and cross-matched (needs ``max_pix_value`` and a
    V column).
    """
    vcol = next((c for c in ("v_jkc_mag", "v_mag") if c in sources.colnames), None)
    if "max_pix_value" not in sources.colnames or vcol is None:
        return float("nan")
    peak = np.asarray(sources["max_pix_value"], dtype=float)
    v = np.asarray(sources[vcol], dtype=float)
    matched = ~np.ma.getmaskarray(sources[vcol])
    sat = matched & np.isfinite(peak) & np.isfinite(v) & (peak >= sat_level)
    return float(np.max(v[sat])) if sat.any() else float("nan")


def sky_surface_brightness(pedestal, zeropoint, pixscale, bias=0.0):
    """Sky surface brightness (mag/arcsec^2) from the background pedestal.

    ``mu = ZP - 2.5 * log10((pedestal - bias) / pixscale^2)`` -- the bias-subtracted
    sky level per pixel, spread to per-arcsec^2 and put through the zero point.
    Returns nan if the sky level is non-positive.

    Caveats: this is a *V-equivalent* brightness, calibrated to Gaia synthetic V
    through the broad IRCUT band, so it carries a sensor- and band-dependent offset
    (compare within a unit, not blindly across S50 vs S30pro); and the ZP folds in
    the ~0.1 mag aperture correction.
    """
    sky = pedestal - bias
    if not np.isfinite(sky) or sky <= 0:
        return float("nan")
    return float(zeropoint - 2.5 * np.log10(sky / pixscale ** 2))


def effective_bortle(sqm):
    """Approximate Bortle class (1-9) from a V-band sky surface brightness."""
    if not np.isfinite(sqm):
        return float("nan")
    for lower, cls in _BORTLE_SQM:
        if sqm >= lower:
            return cls
    return 9
