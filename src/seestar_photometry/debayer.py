"""Bayer demosaic for raw Seestar sub-exposures.

A Seestar *stack* arrives already debayered as three colour planes. A raw sub does not:
it is a single 2-D frame carrying the sensor's colour-filter mosaic, tagged by the
``BAYERPAT`` header keyword (``GRBG`` on the S50 and S30pro).

:func:`debayer` turns that into the ``(3, ny, nx)`` R, G, B cube the rest of the package
expects, so a raw sub and an on-board stack are interchangeable downstream.

The demosaic is deliberately plain bilinear. Anything cleverer (gradient-corrected,
edge-aware) is tuned for photographs of edges, whereas these frames are point sources on
a flat background, where the fancier methods buy nothing measurable and make the noise
harder to reason about.

Pattern orientation, verified on a real S50 sub rather than assumed: for ``GRBG`` as the
array is indexed, green sits on the ``(0,0)`` and ``(1,1)`` pixels of each 2x2 quad, red
on ``(0,1)`` and blue on ``(1,0)``. The check is that the two green sub-lattices agree --
they measured 4 ADU apart with matching sigma, while ``(0,1)`` read 986 and ``(1,0)``
1103. There is no row flip.
"""

import numpy as np

#: Offsets of each colour within the 2x2 mosaic quad, as the array is indexed.
#: ``{pattern: {band: ((row, col), ...)}}``.
BAYER_PATTERNS = {
    "GRBG": {"G": ((0, 0), (1, 1)), "R": ((0, 1),), "B": ((1, 0),)},
    "RGGB": {"R": ((0, 0),), "G": ((0, 1), (1, 0)), "B": ((1, 1),)},
    "BGGR": {"B": ((0, 0),), "G": ((0, 1), (1, 0)), "R": ((1, 1),)},
    "GBRG": {"G": ((0, 0), (1, 1)), "B": ((0, 1),), "R": ((1, 0),)},
}

#: The Seestar's pattern, used when a header omits ``BAYERPAT``.
DEFAULT_PATTERN = "GRBG"

#: 3x3 bilinear kernel. Green is sampled on a quincunx rather than a square grid, but the
#: same kernel serves every plane once the result is divided by the convolved sample mask.
_KERNEL = np.array([[1.0, 2.0, 1.0],
                    [2.0, 4.0, 2.0],
                    [1.0, 2.0, 1.0]], dtype=np.float32)


def pattern_of(header, default=DEFAULT_PATTERN):
    """The ``BAYERPAT`` of a header, upper-cased, falling back to ``default``."""
    value = str(header.get("BAYERPAT", "") or "").strip().upper()
    return value if value in BAYER_PATTERNS else default


def is_bayer(header, data):
    """Whether this HDU looks like an undemosaiced Bayer frame.

    True for 2-D data carrying a ``BAYERPAT``. The keyword matters: a 2-D FITS image
    could be anything, and silently demosaicing a mono frame would quietly triple a
    single channel into a grey "colour" image rather than failing.
    """
    return (data is not None and getattr(data, "ndim", 0) == 2
            and bool(str(header.get("BAYERPAT", "") or "").strip()))


def debayer(raw, pattern=DEFAULT_PATTERN):
    """Bilinear demosaic of a 2-D Bayer frame into a ``(3, ny, nx)`` float32 cube.

    Each colour's samples are scattered onto a full-size grid and convolved, then divided
    by the same convolution of the sample *mask*. Normalising by the mask -- rather than
    using fixed 1/4 and 1/2 weights -- makes the frame edges come out right for free,
    which matters because sources near the edge are still measured.
    """
    from scipy.ndimage import convolve

    raw = np.asarray(raw)
    if raw.ndim != 2:
        raise ValueError(f"expected a 2-D Bayer frame, got shape {raw.shape}")
    if pattern not in BAYER_PATTERNS:
        raise ValueError(f"unknown Bayer pattern {pattern!r}; "
                         f"expected one of {sorted(BAYER_PATTERNS)}")
    offsets = BAYER_PATTERNS[pattern]

    a = raw.astype(np.float32)
    ny, nx = a.shape
    out = np.empty((3, ny, nx), dtype=np.float32)
    for index, band in enumerate(("R", "G", "B")):
        samples = np.zeros((ny, nx), dtype=np.float32)
        mask = np.zeros((ny, nx), dtype=np.float32)
        for (i, j) in offsets[band]:
            samples[i::2, j::2] = a[i::2, j::2]
            mask[i::2, j::2] = 1.0
        num = convolve(samples, _KERNEL, mode="nearest")
        den = convolve(mask, _KERNEL, mode="nearest")
        out[index] = num / np.maximum(den, 1e-6)
    return out


def green_half(raw, pattern=DEFAULT_PATTERN):
    """Half-resolution green plane, by averaging the two green sub-lattices.

    This is the natural image to *register* on: it involves no interpolation, so it
    carries no demosaic artefacts, it is a quarter of the pixels, and green is the
    deepest channel. Used by :mod:`stacking`.
    """
    raw = np.asarray(raw).astype(np.float32)
    greens = BAYER_PATTERNS[pattern]["G"]
    planes = [raw[i::2, j::2] for (i, j) in greens]
    # The two green lattices are the same shape for an even-sized frame; guard the odd
    # case by trimming to the common shape rather than raising.
    ny = min(p.shape[0] for p in planes)
    nx = min(p.shape[1] for p in planes)
    return sum(p[:ny, :nx] for p in planes) / float(len(planes))


def channel_medians(raw, pattern=DEFAULT_PATTERN):
    """Median of each colour's *native* samples, with no interpolation.

    A quick way to see the sensor's channel balance. A raw Seestar sub reads roughly
    R/G/B = 986/965/1103 ADU, while an on-board stack reads ~963/965/964 -- the on-board
    stacker balances the channel backgrounds. Green, which carries the photometry, is
    untouched by that balancing, so it does not affect the calibration either way.
    """
    a = np.asarray(raw).astype(np.float32)
    out = {}
    for band, offsets in BAYER_PATTERNS[pattern].items():
        vals = np.concatenate([a[i::2, j::2].ravel() for (i, j) in offsets])
        out[band] = float(np.median(vals))
    return out
