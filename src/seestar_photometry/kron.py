"""Kron/AUTO elliptical photometry -- a cross-check, not the standard path.

Kept because it is useful for star/galaxy separation and as an independent
total-flux comparison against the fixed-aperture photometry in
:mod:`photometry`. It is deliberately **not** the standard path: its per-source
adaptive aperture introduces source-dependent (and frame-dependent) systematics
that do not cancel cleanly in differential photometry or the zero-point fit. See
``docs/photometry-design.md``.

A known feature: sources whose Kron radius is too small to integrate fall back to a
fixed minimum circular aperture, forming a distinct ``kron_fallback``
sub-population -- visible as a secondary streak in the magnitude-SNR diagram. Filter
on the ``kron_fallback`` column before using Kron fluxes as total magnitudes.
"""

import numpy as np
import sep
from astropy.table import Table, vstack

from .frames import BANDS
from .photometry import SNR_MIN, Extraction, _detect


def _kron_flux(data, objs, err, k=2.5, r_min=1.75, bkgann=None):
    """SExtractor-style Kron (AUTO) elliptical photometry.

    Measures flux within an ellipse of ``k * kron_radius``, falling back to a
    circular aperture of radius ``r_min`` for sources whose Kron radius is too small
    to integrate reliably. Returns ``(flux, fluxerr, flag, fallback)``, where
    ``fallback`` marks the sources measured in the fixed circle.
    """
    x, y, a, b, theta = objs["x"], objs["y"], objs["a"], objs["b"], objs["theta"]
    # SEP's ellipse routines require |theta| strictly < pi/2; extract can return
    # values a float epsilon outside that range, so clip just inside.
    theta = np.clip(theta, -np.pi / 2 + 1e-6, np.pi / 2 - 1e-6)
    kronrad, krflag = sep.kron_radius(data, x, y, a, b, theta, 6.0)
    flux = np.zeros(len(x))
    fluxerr = np.zeros(len(x))
    flag = np.zeros(len(x), dtype=np.int32)
    use_circle = ~(kronrad * np.sqrt(a * b) >= r_min)
    ell = ~use_circle
    if np.any(ell):
        flux[ell], fluxerr[ell], flag[ell] = sep.sum_ellipse(
            data, x[ell], y[ell], a[ell], b[ell], theta[ell],
            k * kronrad[ell], err=err, bkgann=bkgann, subpix=1,
        )
    if np.any(use_circle):
        flux[use_circle], fluxerr[use_circle], flag[use_circle] = sep.sum_circle(
            data, x[use_circle], y[use_circle], r_min,
            err=err, bkgann=bkgann, subpix=1,
        )
    flag |= krflag
    return flux, fluxerr, flag, use_circle


def extract_kron(frame, thresh=2.0, bkgann=None):
    """Kron/AUTO photometry across the R, G, B planes (cross-check only).

    Mirrors :func:`photometry.extract_sources` but with adaptive elliptical
    apertures. The ``sources`` table gains a ``kron_fallback`` boolean column. Only
    sources with SNR > 5 are kept.
    """
    background = np.empty_like(frame.data)
    rms = np.empty(3, dtype=np.float32)
    tables = []
    for i, band in enumerate(BANDS):
        bkg, data_sub, objs = _detect(frame.data[i], thresh)
        background[i] = bkg.back()
        rms[i] = bkg.globalrms
        flux, fluxerr, flag, fallback = _kron_flux(
            data_sub, objs, bkg.globalrms, bkgann=bkgann
        )
        snr = flux / fluxerr
        keep = snr > SNR_MIN
        t = Table()
        t["band"] = np.full(int(keep.sum()), band)
        t["x"] = objs["x"][keep]
        t["y"] = objs["y"][keep]
        t["flux"] = flux[keep]
        t["fluxerr"] = fluxerr[keep]
        t["snr"] = snr[keep]
        t["flag"] = flag[keep]
        t["kron_fallback"] = fallback[keep]
        t["a"] = objs["a"][keep]
        t["b"] = objs["b"][keep]
        t["theta"] = objs["theta"][keep]
        tables.append(t)
    return Extraction(
        background=background, rms=rms, sources=vstack(tables), frame=frame
    )
