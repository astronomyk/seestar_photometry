"""Per-frame quality metrics -- one row of the frame table.

Every frame that survives calibration contributes one row combining its zero point,
photometric scatter, PSF, sky background and source counts. The accumulated table
(``frames.ecsv``) is the thing you query to decide which frames to trust, how deep
the dataset goes, and how conditions varied across a night.

Reading the important columns:

``rms``
    Scatter of the zero-point fit residuals (mag). The master quality number.
    Below ~0.06 is photometric-grade.
``chi2_red``
    Residual scatter against the photon-noise prediction. ~1 would mean pure
    measurement noise; in practice even the best Seestar frames land at 100-200,
    which is the signature of a ~0.03 mag systematic floor, not a broken fit. Use it
    to spot frames that are *much* worse than their peers, not as an absolute test.
``v_lim_5sigma`` / ``v_lim_100sigma``
    Faint and bright-ish ends of the usable range, extrapolated from the zero point
    and the aperture noise. Optimistic in absolute terms -- see
    :func:`calibration.limiting_mag`.
``v_sat``
    Empirical saturation limit: brighter than this the flux stops tracking V.
``sky_sb`` / ``bortle``
    V-equivalent sky brightness and its Bortle class. Comparable within a unit;
    across models it carries a sensor-dependent offset.

See ``docs/frame-table.md``.
"""

from pathlib import Path

import numpy as np

from . import calibration


def frame_quality(extraction, cal, provenance=None):
    """One row of per-frame quality metrics, as a flat dict.

    Parameters
    ----------
    extraction : Extraction
        A cross-matched extraction (``match_gaia`` already applied).
    cal : Calibration
        The green-band zero-point fit for the same frame.
    provenance : callable, optional
        ``provenance(frame) -> dict`` of extra columns to merge in -- dataset name,
        stacking manifest fields, binning, whatever a particular project tracks.
        Keeps project-specific bookkeeping out of the core schema without needing a
        second table builder.
    """
    from . import astrometry, frames, photometry

    frame = extraction.frame
    meta = frames.frame_metadata(frame)
    g = extraction.band("G")

    vcol = next((c for c in ("v_jkc_mag", "v_mag") if c in g.colnames), None)
    n_matched = int((~np.ma.getmaskarray(g[vcol])).sum()) if vcol else 0
    # Aperture noise is the (constant) green aperture flux error; the median is
    # robust to the few apertures truncated at the frame edge.
    sigma_aper = float(np.median(np.asarray(g["fluxerr"], dtype=float)))

    bg = photometry.fit_background(extraction, band="G")
    pixscale = astrometry.pixel_scale(frame)
    sky_sb = calibration.sky_surface_brightness(
        bg["pedestal"], cal.zeropoint, pixscale, frame.header.get("BIAS", 0.0)
    )

    row = {
        # identity / provenance
        "path": str(frame.path),
        "frame": Path(frame.path).name,
        "layout": frame.layout,
        "model": frame.model,
        "unit": meta["unit"],
        "telescope": str(frame.header.get("TELESCOP", "?")),
        "object": str(frame.header.get("OBJECT", "")),
        "filter": str(frame.header.get("FILTER", "")),
        # timing / geometry
        "date_obs": str(frame.header.get("DATE-OBS", "")),
        "obs_start": meta["obs_start"] or "",
        "obs_end": meta["obs_end"] or "",
        "eqmode": meta["eqmode"],  # 0 = Alt-Az, 1 = equatorial
        "n_exp": meta["n_exp"],
        "exptime": meta["exptime"],
        "total_exptime": meta["total_exptime"],
        "airmass": meta["airmass"],
        "site_lat": meta["site_lat"],
        "site_lon": meta["site_lon"],
        "ccd_temp": meta["ccd_temp"],
        "pixscale": pixscale,
        # counts
        "n_sources": len(extraction.sources),
        "n_green": len(g),
        "n_matched": n_matched,
        "n_cal": cal.n_stars,
        # calibration
        "zeropoint": cal.zeropoint,
        "zeropoint_err": cal.zeropoint_err,
        "colour_term": cal.colour_term,
        "colour_term_err": cal.colour_term_err,
        "colour0": cal.colour0,
        "rms": cal.rms,
        "chi2_red": cal.chi2_red,
        "fit_mag_lo": cal.mag_range[0] if cal.mag_range else np.nan,
        "fit_mag_hi": cal.mag_range[1] if cal.mag_range else np.nan,
        # depth / dynamic range
        "v_lim_5sigma": calibration.limiting_mag(cal.zeropoint, sigma_aper),
        "v_lim_100sigma": calibration.limiting_mag(cal.zeropoint, sigma_aper, nsigma=100.0),
        "v_sat": calibration.saturation_mag(g),
        "sigma_aper": sigma_aper,
        # PSF
        "fwhm_R": float(extraction.fwhm[0]),
        "fwhm_G": float(extraction.fwhm[1]),
        "fwhm_B": float(extraction.fwhm[2]),
        "aperture_R": float(extraction.aperture[0]),
        "aperture_G": float(extraction.aperture[1]),
        "aperture_B": float(extraction.aperture[2]),
        # sky
        "sky_R": float(extraction.rms[0]),
        "sky_G": float(extraction.rms[1]),
        "sky_B": float(extraction.rms[2]),
        "sky_pedestal": bg["pedestal"],
        "bg_poly": bg["coeffs"],
        "bg_resid_std": bg["resid_std"],
        "sky_sb": sky_sb,
        "bortle": calibration.effective_bortle(sky_sb),
    }
    if provenance is not None:
        row.update(provenance(frame))
    return row


def onboard_quality(frame):
    """The server's own quality metrics, lifted from a CrowdSky header.

    CrowdSky runs an earlier generation of this same pipeline server-side and records
    the result in the primary header (``ZPG``, ``ZPSCTG``, ``FWHMG``, ``SKYRMSG``,
    ``SQMPHOT``, ...). Those numbers are *not* used for science here -- we re-measure
    everything -- but comparing them against our own row is a cheap, independent
    check that a frame was reduced sensibly. Returns ``{}`` for a native frame.
    """
    h = frame.header
    keys = {
        "onboard_zp_R": "ZPR", "onboard_zp_G": "ZPG", "onboard_zp_B": "ZPB",
        "onboard_rms_R": "ZPSCTR", "onboard_rms_G": "ZPSCTG", "onboard_rms_B": "ZPSCTB",
        "onboard_ct_R": "ZPCTR", "onboard_ct_G": "ZPCTG", "onboard_ct_B": "ZPCTB",
        "onboard_fwhm_R": "FWHMR", "onboard_fwhm_G": "FWHMG", "onboard_fwhm_B": "FWHMB",
        "onboard_sky_G": "SKYRMSG", "onboard_sqm": "SQMPHOT",
        "onboard_n_cal": "ZPNSTAR", "onboard_pixscale": "PIXSCALE",
    }
    out = {k: h[v] for k, v in keys.items() if v in h}
    return out


def photometric(frames, rms_max=0.06):
    """Subset a frame table to photometric-grade rows (``rms < rms_max``)."""
    return frames[np.asarray(frames["rms"], dtype=float) < rms_max]
