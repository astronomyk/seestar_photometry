"""Extended-emission contamination for targets sitting on a host.

A supernova in a nearby galaxy, or a star projected on a nebula, has host light
inside its aperture. Fixed-aperture photometry attributes that light to the target,
so the measured magnitude is too bright -- badly so on a bright host, and by an amount
that changes with seeing and pointing, which is worse than a constant offset because
it adds structure to the light curve.

The approach here is **empirical, not model-based**: sample the aperture flux at the
target's galactocentric radius over clean azimuths and take the median as the host
contribution. That works where a smooth (Sersic or similar) model fails, because real
hosts have structure at a given radius -- NGC 3310's circumnuclear starburst ring sits
almost exactly at the radius of SN 2026sqf, and a smooth model under-predicts it
substantially.

Assumptions worth checking before trusting the result: the host is roughly
axisymmetric at the target's radius, and enough azimuths are clean of stars. Both are
reported back so you can tell when they fail.
"""

import numpy as np


def star_mask(data_sub, rms, nucleus_xy, thresh=5.0, keep_radius=6.0, grow=3.5):
    """Mask of field stars, deliberately excluding the host's own core.

    Detections within ``keep_radius`` pixels of ``nucleus_xy`` are *not* masked: the
    galaxy nucleus is itself a detection, and masking it (plus its extended wings)
    would remove exactly the emission we are trying to measure.

    ``grow`` scales each detection's ellipse, so bright stars are masked out to their
    wings rather than just their cores.
    """
    import sep

    objs = sep.extract(data_sub, thresh, err=rms)
    xn, yn = nucleus_xy
    dist = np.hypot(objs["x"] - xn, objs["y"] - yn)
    mask = np.zeros(data_sub.shape, bool)
    for i, o in enumerate(objs):
        if dist[i] >= keep_radius:
            sep.mask_ellipse(
                mask, o["x"], o["y"], o["a"], o["b"], o["theta"], r=grow
            )
    return mask


def galaxy_contamination(data_sub, nucleus_xy, target_xy, aperture, mask=None,
                         step_deg=15.0, min_azimuths=4):
    """Empirical host flux inside the target aperture (ADU).

    Places the same aperture at the target's galactocentric radius around the nucleus,
    skipping azimuths within ``2 * aperture`` of the target itself and any that land on
    a masked star, and takes the median.

    Returns a dict with ``adu`` (the median host flux -- subtract this from the target
    flux), ``std`` (azimuthal scatter, a fair estimate of the systematic uncertainty
    on that subtraction), ``radius_px``, ``n_azimuth`` (how many clean samples
    contributed) and ``samples``. ``adu`` is 0.0 with ``n_azimuth`` below
    ``min_azimuths``, so a hopelessly crowded field degrades to no correction rather
    than to a wild one -- check ``n_azimuth`` before trusting the number.
    """
    import sep

    xn, yn = float(nucleus_xy[0]), float(nucleus_xy[1])
    xs, ys = float(target_xy[0]), float(target_xy[1])
    radius = np.hypot(xs - xn, ys - yn)
    ny, nx = data_sub.shape
    edge = int(np.ceil(aperture)) + 1

    samples = []
    for angle in np.deg2rad(np.arange(0.0, 360.0, step_deg)):
        px, py = xn + radius * np.cos(angle), yn + radius * np.sin(angle)
        if np.hypot(px - xs, py - ys) < 2 * aperture:
            continue  # too close to the target itself
        iy, ix = int(py), int(px)
        if not (edge <= iy < ny - edge and edge <= ix < nx - edge):
            continue  # aperture would run off the frame
        if mask is not None and mask[iy - edge:iy + edge + 1,
                                    ix - edge:ix + edge + 1].any():
            continue  # a field star sits here
        flux, _, _ = sep.sum_circle(
            data_sub, np.array([px]), np.array([py]), aperture
        )
        samples.append(float(flux[0]))

    samples = np.asarray(samples, dtype=float)
    enough = len(samples) >= min_azimuths
    return {
        "adu": float(np.median(samples)) if enough else 0.0,
        "std": float(samples.std()) if len(samples) > 1 else float("nan"),
        "radius_px": float(radius),
        "n_azimuth": int(len(samples)),
        "samples": samples,
    }


def cutout(data_sub, nucleus_xy, target_xy, aperture, half=70):
    """A postage stamp around host and target, for the diagnostic figure.

    Returns the ``stamp`` plus ``target_xy``/``nucleus_xy`` in stamp coordinates, the
    ``aperture`` and the galactocentric ``radius_px``, so the figure can draw the
    aperture and the azimuthal sampling ring that produced the contamination estimate.
    """
    xn, yn = float(nucleus_xy[0]), float(nucleus_xy[1])
    xs, ys = float(target_xy[0]), float(target_xy[1])
    ny, nx = data_sub.shape
    y0, y1 = max(int(yn) - half, 0), min(int(yn) + half, ny)
    x0, x1 = max(int(xn) - half, 0), min(int(xn) + half, nx)
    return {
        "stamp": data_sub[y0:y1, x0:x1].astype(np.float32),
        "target_xy": np.array([xs - x0, ys - y0], dtype=float),
        "nucleus_xy": np.array([xn - x0, yn - y0], dtype=float),
        "aperture": float(aperture),
        "radius_px": float(np.hypot(xs - xn, ys - yn)),
    }
