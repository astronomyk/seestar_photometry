"""System characterisation: limiting magnitude vs exposure time.

Each frame is one stack of ``n_exp`` sub-exposures, so its 5-sigma limit (the
``v_lim_5sigma`` column of the frame table) corresponds to that frame's
``total_exptime`` -- a non-round value (e.g. 390 or 560 s), not a tidy 15 minutes.
To compare frames or quote a headline depth they must first be put on a common
exposure.

For a background-limited measurement the 5-sigma flux scales as 1/sqrt(t), so the
limiting magnitude deepens by ``1.25 * log10(t2/t1)``. This module rescales per-frame
limits to a reference exposure and builds the limit-vs-exposure relation per group
(instrument model, individual unit, dataset, ...).

Caveats, both of which matter before quoting a number:

* The sqrt(t) law is exact only while the co-add stays background-limited.
  Extrapolating to 1-16 h assumes no systematic floor (flat-field residuals,
  imperfect alignment). This module applies the *theoretical* scaling; it does not
  co-add frames, so long-exposure values should be validated by actually stacking.
* ``v_lim_5sigma`` is itself a noise-model extrapolation of the bright-star zero
  point (see :func:`calibration.limiting_mag`), and empirically runs optimistic by
  roughly 0.7-1.1 mag against what is actually recoverable in a stack. Treat it as
  a relative quality metric, not an absolute detection promise.
"""

import numpy as np
from astropy.table import Table

#: Reference exposure for the headline limit: 15 minutes, in seconds.
EXPTIME_REF = 900.0

#: Standard exposure grid for the limit-vs-time relation (s): 1 min .. 16 h.
EXPTIME_GRID = (60.0, 900.0, 3600.0, 14400.0, 57600.0)

#: Calibration scatter (mag) below which a frame counts as photometric-grade.
RMS_PHOTOMETRIC = 0.06


def scale_limit(v_lim, t_from, t_to=EXPTIME_REF):
    """Scale a 5-sigma limiting magnitude from exposure ``t_from`` to ``t_to``.

    Background-limited: the 5-sigma flux goes as 1/sqrt(t), so the limit deepens by
    ``1.25 * log10(t_to / t_from)``. Arrays broadcast.
    """
    return np.asarray(v_lim, float) + 1.25 * np.log10(
        np.asarray(t_to, float) / np.asarray(t_from, float)
    )


def detection_limit(frames, t_ref=EXPTIME_REF, rms_max=RMS_PHOTOMETRIC, group="model"):
    """Per-group 5-sigma detection limit at a single reference exposure.

    Each photometric-grade frame's ``v_lim_5sigma`` (at its own ``total_exptime``) is
    rescaled to ``t_ref`` and summarised by ``group``. Returns a table with the
    median limit, its scatter and the frame count.
    """
    q = frames[np.asarray(frames["rms"], float) < rms_max]
    v_ref = scale_limit(q["v_lim_5sigma"], q["total_exptime"], t_ref)
    keys = np.asarray(q[group])
    rows = []
    for key in sorted(set(keys.tolist())):
        sel = keys == key
        rows.append({
            group: key,
            "exptime": float(t_ref),
            "n_frames": int(sel.sum()),
            "v_lim_median": float(np.median(v_ref[sel])),
            "v_lim_std": float(np.std(v_ref[sel])),
        })
    return Table(rows)


def limit_vs_exptime(frames, exptimes=EXPTIME_GRID, rms_max=RMS_PHOTOMETRIC,
                     group="model"):
    """Expected 5-sigma limit across an exposure grid, per group.

    For each exposure every photometric-grade frame's limit is scaled from its own
    ``total_exptime``, then summarised by ``group``. Traces limiting magnitude vs
    integration time (slope 1.25 per dex, by construction). Returns a long-format
    table ``[group, exptime, v_lim_median, v_lim_std, n_frames]``.

    See the module caveats before quoting the long-exposure end.
    """
    q = frames[np.asarray(frames["rms"], float) < rms_max]
    keys = np.asarray(q[group])
    t_frame = np.asarray(q["total_exptime"], float)
    v = np.asarray(q["v_lim_5sigma"], float)
    rows = []
    for key in sorted(set(keys.tolist())):
        sel = keys == key
        for t in exptimes:
            v_t = scale_limit(v[sel], t_frame[sel], t)
            rows.append({
                group: key,
                "exptime": float(t),
                "v_lim_median": float(np.median(v_t)),
                "v_lim_std": float(np.std(v_t)),
                "n_frames": int(sel.sum()),
            })
    return Table(rows)


def fit_depth_model(frames, rms_max=RMS_PHOTOMETRIC):
    """Fit ``v_lim = a + b*log10(t) + c*sky_sb + d*fwhm`` to a frame table.

    A *noise* model, not a zero-point correction: it describes how deep a frame goes
    given its integration, sky brightness and seeing. Use it to predict the depth of
    a planned exposure, or to compare units after removing the conditions they
    happened to observe in.

    The coefficients are per-unit in practice -- each physical Seestar responds
    differently to sky and seeing -- so fit each unit separately, and be wary of a
    unit whose data span only a narrow range of conditions: there ``c`` and ``d`` are
    poorly constrained and will happily take on implausible values.

    Returns a dict with ``coeffs`` (a, b, c, d), ``rms`` of the residuals and
    ``n_frames``.
    """
    q = frames[np.asarray(frames["rms"], float) < rms_max]
    t = np.asarray(q["total_exptime"], float)
    sky = np.asarray(q["sky_sb"], float)
    fwhm = np.asarray(q["fwhm_G"], float)
    v = np.asarray(q["v_lim_5sigma"], float)
    ok = np.isfinite(t) & (t > 0) & np.isfinite(sky) & np.isfinite(fwhm) & np.isfinite(v)
    design = np.column_stack([
        np.ones(int(ok.sum())), np.log10(t[ok]), sky[ok], fwhm[ok],
    ])
    coeffs, *_ = np.linalg.lstsq(design, v[ok], rcond=None)
    resid = v[ok] - design @ coeffs
    return {
        "coeffs": coeffs,
        "rms": float(np.std(resid)),
        "n_frames": int(ok.sum()),
    }
