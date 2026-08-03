"""Light-curve construction: timing, comparison selection, ensemble differential.

The data model is two tables, built by :mod:`pipeline`:

``stars.ecsv``
    One row per catalogue source in the field -- position, magnitudes, colour,
    variability flag, separation from the target. The small table you *query* to pick
    a comparison ensemble.

``measurements.ecsv``
    The long/tidy table, one row per ``(source_id, frame, band)``: forced-aperture
    flux, instrumental and calibrated magnitude, mid-exposure time and airmass.

Two choices make the <1% floor reachable and are worth understanding before
changing anything here:

**Forced photometry, not detect-then-match.** Every catalogue position is measured in
every frame whether or not a source was detected there. A star therefore keeps a row
in every frame it lands on, so series never go ragged and the comparison ensemble is
*identical* frame to frame. With detect-then-match, a faint comparison drops out
exactly on the frames where conditions were poor -- precisely when you need it -- and
the ensemble composition changes underneath the light curve.

**Each comparison referenced to its own catalogue magnitude.** The per-frame zero
point is the mean over comparisons of ``(catalogue mag - instrumental mag)``, not a
mean flux ratio. That way a comparison dropping in or out between frames does not
shift the zero point, and the scatter of the per-comparison zero points is the
genuine measurement error rather than the intrinsic brightness spread of the
ensemble.

See ``docs/light-curves.md``.
"""

import warnings

import numpy as np


def frame_times(header, target, total_exptime=None):
    """Mid-exposure timestamps for a stacked frame.

    Getting this right matters as much as the photometry for anything with structure
    on tens of minutes: a stack is a 10-15 minute integration, so a mid-point error
    of half that smears and shifts a transit or a short-period pulsation.

    Two cases, in order of preference:

    1. The header records the true exposure span (``OB-START``/``OB-END``, written by
       CrowdSky). The mid-point is then exact, including inter-sub overhead.
    2. Only ``DATE-OBS`` and the on-sky integration are known (native Seestar). The
       mid-point is taken as ``DATE-OBS + total_exptime/2``. This ignores overhead
       between subs -- unknowable from the header -- so it runs slightly early; for a
       390 s on-sky stack spanning ~650 s of wall clock the bias is ~2 minutes.

    The barycentric correction uses the target position: over a sub-degree field the
    per-star light-travel-time difference is < 1 s, so one coordinate serves the whole
    frame.

    Parameters
    ----------
    header : fits.Header
        Needs ``DATE-OBS`` (or ``OB-START``); ``SITELAT``/``SITELONG`` for BJD.
    target : astropy.coordinates.SkyCoord
        Target (field) coordinate, for the barycentric correction.
    total_exptime : float, optional
        On-sky integration (s), used only in case 2.

    Returns
    -------
    dict
        ``mjd_obs`` (UTC MJD at the start epoch), ``mjd_mid`` (UTC MJD at mid
        exposure), ``bjd_tdb`` (barycentric Julian date, TDB, at mid exposure; nan if
        the site is unknown) and ``time_source`` (``"span"`` or ``"exptime"``).
    """
    import astropy.units as u
    from astropy.coordinates import EarthLocation
    from astropy.time import Time
    from astropy.utils import iers

    iers.conf.auto_download = False
    # Tolerate a stale/extrapolated bundled table (offline, or observation dates
    # beyond what the shipped table covers).
    iers.conf.auto_max_age = None

    from .frames import exposure_span

    start, end = exposure_span(header)
    t_obs = Time(start, format="isot", scale="utc")
    if end is not None:
        t_end = Time(end, format="isot", scale="utc")
        t_mid = t_obs + (t_end - t_obs) / 2
        source = "span"
    else:
        half = 0.0 if (total_exptime is None or not np.isfinite(total_exptime)) \
            else float(total_exptime) / 2.0
        t_mid = t_obs + half * u.s
        source = "exptime"

    lat, lon = header.get("SITELAT"), header.get("SITELONG")
    bjd_tdb = float("nan")
    if lat is not None and lon is not None:
        loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ltt = t_mid.light_travel_time(target, "barycentric", location=loc)
            bjd_tdb = float((t_mid.tdb + ltt).jd)
    return {
        "mjd_obs": float(t_obs.mjd),
        "mjd_mid": float(t_mid.mjd),
        "bjd_tdb": bjd_tdb,
        "time_source": source,
    }


#: Catalogue columns carried into the stars table.
STAR_COLUMNS = (
    "source_id", "ra", "dec",
    "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp",
    "v_jkc_mag", "b_jkc_mag", "r_jkc_mag", "phot_variable_flag",
)


def build_stars(catalogue, source_ids, target_radec, target_id=None):
    """One row per catalogue source that was measured, with separation from the target.

    The nearest source to ``target_radec`` is flagged ``is_target`` unless an explicit
    ``target_id`` is given. Prefer passing ``target_id`` when you know it: "nearest to
    the pointing" is right for a field centred on its target but wrong for a target
    that is not the brightest thing at the centre, or a crowded field.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.table import Table

    present = np.isin(np.asarray(catalogue["source_id"]), np.asarray(source_ids))
    sub = catalogue[present]
    stars = Table()
    for col in STAR_COLUMNS:
        if col in sub.colnames:
            stars[col] = sub[col]
    target = SkyCoord(target_radec[0] * u.deg, target_radec[1] * u.deg)
    coords = SkyCoord(np.asarray(sub["ra"]) * u.deg, np.asarray(sub["dec"]) * u.deg)
    stars["sep_target_arcmin"] = coords.separation(target).arcmin
    stars["is_target"] = np.zeros(len(stars), dtype=bool)
    if len(stars):
        if target_id is not None:
            hit = np.asarray(stars["source_id"]) == target_id
            if not hit.any():
                raise ValueError(f"target {target_id} was not measured in any frame")
            stars["is_target"][np.argmax(hit)] = True
        else:
            stars["is_target"][int(np.argmin(stars["sep_target_arcmin"]))] = True
    stars.sort("sep_target_arcmin")
    return stars


def _target_row(stars, target_id):
    """The stars-table row for ``target_id`` (or the ``is_target`` flag)."""
    if target_id is None:
        if "is_target" not in stars.colnames or not stars["is_target"].any():
            raise ValueError("no target_id given and no is_target flag set")
        return stars[np.asarray(stars["is_target"])][0]
    match = stars[np.asarray(stars["source_id"]) == target_id]
    if not len(match):
        raise ValueError(f"target {target_id} not in stars table")
    return match[0]


def target_id_of(stars):
    """The ``source_id`` of the flagged target."""
    return int(_target_row(stars, None)["source_id"])


def select_comparisons(
    stars, target_id=None, max_sep_arcmin=None, dmag=1.0, mag_range=None,
    colour_tol=0.3, exclude_variable=True, mag_col="v_jkc_mag", colour_col="bp_rp",
):
    """Choose comparison-star ``source_id``\\ s from the stars table.

    Filters on what actually matters for differential photometry:

    * **proximity** -- shared atmosphere, and on an Alt-Az mount a shared amount of
      field rotation, so the position-dependent aperture loss largely cancels;
    * **similar brightness** -- similar noise, and away from both saturation and the
      detection floor;
    * **similar colour** -- minimises differential-extinction and colour-term
      residuals;

    and drops catalogue-flagged variables. All cuts are optional; pass ``None`` to
    skip one.

    A larger ensemble is not automatically better: adding faint or distant stars adds
    noise and systematics faster than it averages them down. Somewhere around 10-30
    well-matched comparisons is usually the sweet spot, and the per-comparison
    diagnostic curves from :func:`comparison_curves` are how you check the ones you
    picked are behaving.

    Returns the selected rows (the target is always excluded).
    """
    tgt = _target_row(stars, target_id)
    keep = np.asarray(stars["source_id"]) != tgt["source_id"]

    mag = np.asarray(stars[mag_col], dtype=float)
    keep &= np.isfinite(mag)
    if mag_range is not None:
        keep &= (mag >= min(mag_range)) & (mag <= max(mag_range))
    elif dmag is not None:
        keep &= np.abs(mag - float(tgt[mag_col])) <= dmag
    if colour_tol is not None and colour_col in stars.colnames:
        col = np.asarray(stars[colour_col], dtype=float)
        keep &= np.isfinite(col) & (np.abs(col - float(tgt[colour_col])) <= colour_tol)
    if max_sep_arcmin is not None and "sep_target_arcmin" in stars.colnames:
        keep &= np.asarray(stars["sep_target_arcmin"], dtype=float) <= max_sep_arcmin
    if exclude_variable and "phot_variable_flag" in stars.colnames:
        flag = np.ma.getdata(stars["phot_variable_flag"]).astype(str)
        keep &= flag != "VARIABLE"
    return stars[keep]


def differential_lightcurve(
    measurements, target_id, comps, band="G", time_col="bjd_tdb",
    inst_col="m_inst", cat_mag_col="v_jkc_mag", min_comp=3, max_flag=0,
):
    """Ensemble differential light curve of one target in one band.

    Each frame gets an **ensemble zero point** -- the mean, over its valid comparison
    stars, of ``(catalogue mag - instrumental mag)``. The target's calibrated
    magnitude in that frame is ``m_target + ensemble ZP``. This removes the per-frame
    transparency (common to target and comparisons, so it cancels) *and* ties the
    target to the comparisons' catalogue scale in one step.

    Only valid measurements are used: on chip, finite magnitude, SEP ``flag <=
    max_flag``.

    Parameters
    ----------
    measurements : astropy.table.Table
        The ``measurements.ecsv`` table.
    target_id : int
        Target ``source_id``.
    comps : astropy.table.Table
        Comparison stars from :func:`select_comparisons`; needs ``source_id`` and the
        ``cat_mag_col`` catalogue magnitude.
    band : str
        Band to build the curve in (``"G"`` ~ Johnson V).
    time_col, inst_col, cat_mag_col : str
        Time, instrumental-magnitude and catalogue-magnitude columns.
    min_comp : int
        Minimum valid comparisons required to keep a frame. Below ~3 the ensemble
        zero point is not averaging anything and its error is unmeasurable.
    max_flag : int
        Maximum SEP flag accepted (0 = clean apertures only).

    Returns
    -------
    astropy.table.Table
        One row per frame: ``time``, ``mag`` (calibrated target magnitude),
        ``mag_err``, ``dmag`` (``mag`` minus its median), ``zp_ens``, ``zp_scatter``,
        ``n_comp``, ``airmass``, ``frame``. Sorted by time; ``meta`` records the
        target, band and ensemble size.
    """
    from astropy.table import Table

    cat = {int(s): float(m_) for s, m_ in zip(comps["source_id"], comps[cat_mag_col])
           if np.isfinite(float(m_))}
    m = measurements[np.asarray(measurements["band"]) == band]
    sid = np.asarray(m["source_id"])
    inst = np.asarray(m[inst_col], dtype=float)
    valid = (
        np.asarray(m["on_chip"])
        & np.isfinite(inst)
        & (np.asarray(m["flag"]) <= max_flag)
    )
    is_comp = np.array([int(s) in cat for s in sid])
    comp_cat = np.array([cat.get(int(s), np.nan) for s in sid])
    frames_col = np.asarray(m["frame"])

    # ``dmag`` is in the schema from the start so an empty curve still has the full set
    # of columns -- otherwise a downstream call fails with a bare KeyError instead of
    # the obvious "no epochs survived", which is the actual problem.
    names = ("time", "mag", "mag_err", "zp_ens", "zp_scatter", "n_comp",
             "airmass", "frame", "dmag")
    rows = []
    for frame in np.unique(frames_col):
        in_frame = frames_col == frame
        t_sel = in_frame & valid & (sid == int(target_id))
        c_sel = in_frame & valid & is_comp
        if not t_sel.any() or int(c_sel.sum()) < min_comp:
            continue
        zp = comp_cat[c_sel] - inst[c_sel]  # per-comparison zero point
        ti = np.where(t_sel)[0][0]
        tgt_err = float(m["mag_err"][ti]) if "mag_err" in m.colnames else 0.0
        # Error on the *mean* zero point, not the spread of the ensemble: the spread
        # is dominated by how well each comparison's catalogue magnitude is known.
        zp_scatter = float(np.std(zp, ddof=1)) if len(zp) > 1 else 0.0
        zp_err = zp_scatter / np.sqrt(len(zp)) if len(zp) > 1 else 0.0
        rows.append((
            float(m[time_col][ti]),
            float(inst[ti] + np.mean(zp)),
            float(np.hypot(tgt_err, zp_err)),
            float(np.mean(zp)),
            zp_scatter,
            int(c_sel.sum()),
            float(m["airmass"][ti]) if "airmass" in m.colnames else np.nan,
            str(frame),
        ))

    if rows:
        # dmag is filled after sorting, so build without it then append.
        lc = Table(rows=[r + (0.0,) for r in rows], names=names)
        lc.sort("time")
        lc["dmag"] = np.asarray(lc["mag"]) - np.median(np.asarray(lc["mag"]))
    else:
        lc = Table(names=names)
    lc.meta.update({
        "target_id": int(target_id), "band": band, "n_comp_pool": len(cat),
        "scatter": float(np.std(np.asarray(lc["mag"]))) if len(lc) else np.nan,
    })
    return lc


def comparison_curves(measurements, comps, band="G", min_comp=3, **kwargs):
    """A differential light curve for each comparison, treating it as the target.

    The single most informative check on an ensemble. Each comparison is measured
    against all the *others*, so a genuinely constant star gives a flat curve at the
    noise floor, while a variable, a blend, a saturated star or one falling off the
    chip stands out immediately -- and would otherwise quietly inflate the target's
    scatter or, worse, imprint its own variability on the target.

    Returns ``{source_id: lightcurve_table}``, each carrying ``meta["scatter"]``.
    Sorting the result by that scatter is the fastest way to find a bad comparison.
    """
    out = {}
    for sid in np.asarray(comps["source_id"]):
        others = comps[np.asarray(comps["source_id"]) != sid]
        if len(others) < min_comp:
            continue
        out[int(sid)] = differential_lightcurve(
            measurements, int(sid), others, band=band, min_comp=min_comp, **kwargs
        )
    return out


def periodogram(lc, min_period=0.01, max_period=None, samples_per_peak=10,
                mag_col="dmag"):
    """Lomb-Scargle periodogram of a light curve.

    Returns ``{"frequency", "power", "best_period", "best_power", "fap"}`` with
    periods in the same time unit as the curve (days for ``bjd_tdb``).

    Watch for aliases: a nightly-cadence dataset has strong power at 1 day and its
    harmonics, and a period suspiciously close to an integer fraction of a day
    deserves the phase-folded plot before it is believed.
    """
    from astropy.timeseries import LombScargle

    t = np.asarray(lc["time"], dtype=float)
    y = np.asarray(lc[mag_col], dtype=float)
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if len(t) < 5:
        raise ValueError(f"need at least 5 finite points for a periodogram, got {len(t)}")
    baseline = t.max() - t.min()
    max_period = max_period if max_period is not None else baseline
    ls = LombScargle(t, y)
    frequency, power = ls.autopower(
        minimum_frequency=1.0 / max_period,
        maximum_frequency=1.0 / min_period,
        samples_per_peak=samples_per_peak,
    )
    best = int(np.argmax(power))
    return {
        "frequency": frequency,
        "power": power,
        "best_period": float(1.0 / frequency[best]),
        "best_power": float(power[best]),
        "fap": float(ls.false_alarm_probability(power[best])),
    }


def phase_fold(lc, period, t0=None, mag_col="dmag"):
    """Phase-fold a light curve. Returns ``(phase, mag)`` with phase in [0, 1)."""
    t = np.asarray(lc["time"], dtype=float)
    t0 = float(np.nanmin(t)) if t0 is None else float(t0)
    return ((t - t0) / period) % 1.0, np.asarray(lc[mag_col], dtype=float)
