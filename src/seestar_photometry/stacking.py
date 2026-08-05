"""Register and co-add raw Seestar sub-exposures into a stack.

The Seestar stacks on board, but only into its own bins, and only if you let every sub
make the round trip to the scope and back. Stacking locally from the raws lets you choose
the bin boundaries, keep every sub, and stack the same photons more than one way.

The pipeline is deliberately plain, and the omissions are the point -- if a local stack
matches an on-board one, the on-board round trip bought nothing:

1. **Bilinear demosaic to full resolution** (:mod:`debayer`), so the output has the same
   shape and pixel scale as an on-board stack and the comparison isolates the *stacking*
   rather than a resampling difference.
2. **Register each sub to the first** with astroalign, as a similarity transform (shift,
   rotation, scale). Rotation is essential, not optional: these are Alt-Az frames and the
   field rotates by several degrees within a 15-minute bin.
3. **Weighted mean, weighted by the warped footprint.** Rotating a frame leaves its
   corners uncovered; without the weight map those corners average against zeros and drag
   the sky level down.

Deliberately *not* done, so any shortfall is attributable: no outlier rejection, no
per-frame quality weighting, no channel balancing, no gradient removal.

Registration is fitted on the half-resolution green plane -- no interpolation, a quarter
of the pixels, and the deepest channel -- then lifted to full resolution.

Needs the ``stack`` extra (``astroalign``, ``scikit-image``).
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io import fits

from .debayer import debayer, green_half, pattern_of

#: astroalign source-detection threshold, in sigma.
DETECTION_SIGMA = 5.0

#: Below this many matched sources the transform is not trustworthy.
MIN_MATCH = 6

#: Median residual (px) above which the transform is not trustworthy.
MAX_RESID_PX = 1.0

#: Fraction of a sub's weight below which an output pixel is left empty.
MIN_WEIGHT = 0.5


@dataclass
class StackReport:
    """What happened during a stack -- read this before trusting the result.

    ``n_bad`` counts subs whose registration was rejected; ``reasons`` says why. A
    ``rot_span`` of several degrees is normal and expected for Alt-Az data over a 15-minute
    bin. ``cover_frac`` well below 1 means field rotation left a large uncovered border,
    which is also normal -- those pixels are zero and excluded by the weight map.
    """

    n_in: int
    n_ok: int
    n_bad: int
    resid_median_px: float
    rot_span_deg: float
    cover_frac: float
    total_exptime: float
    reasons: list = field(default_factory=list)

    def __str__(self):
        return (f"{self.n_ok}/{self.n_in} subs stacked "
                f"({self.total_exptime:.0f} s on sky), "
                f"residual {self.resid_median_px:.2f} px, "
                f"field rotation {self.rot_span_deg:.2f} deg, "
                f"coverage {self.cover_frac:.1%}")


def _scale_transform(tf, factor=2.0):
    """Lift a similarity transform fitted at 1/factor resolution to full resolution.

    Rotation and scale are resolution-independent; only the translation scales. The
    half-res green pixel ``(i, j)`` samples full-res ``(2i+0.5, 2j+0.5)``, and conjugating
    the transform by that mapping leaves a residual offset of ``0.5*(I - R)`` -- about
    1e-3 px for under a degree of rotation, so it is dropped.
    """
    from skimage.transform import SimilarityTransform

    return SimilarityTransform(scale=tf.scale, rotation=tf.rotation,
                               translation=np.asarray(tf.translation) * factor)


def coadd(paths, pattern=None, detection_sigma=DETECTION_SIGMA,
          min_match=MIN_MATCH, max_resid_px=MAX_RESID_PX):
    """Register and co-add raw Bayer subs. Returns ``(cube, weight, header, report)``.

    Streams one sub at a time into a single accumulator, so memory stays flat however many
    subs are passed -- a full night is thousands of frames.

    Parameters
    ----------
    paths : sequence of path-like
        Raw sub-exposures. The **first** is the registration reference and supplies the
        output header, so pass them in a deliberate order (time order is usual).
    pattern : str, optional
        Bayer pattern; read from the first sub's ``BAYERPAT`` if not given.
    detection_sigma, min_match, max_resid_px
        Registration tuning. A sub failing ``min_match`` or ``max_resid_px`` is skipped
        and recorded in the report rather than silently included -- a bad transform smears
        every star in the stack, so rejecting is much better than averaging it in.

    Returns
    -------
    cube : np.ndarray
        ``(3, ny, nx)`` float32 weighted-mean stack, zero where nothing was covered.
    weight : np.ndarray
        ``(ny, nx)`` float64 coverage map: how many subs contributed to each pixel.
    header : fits.Header
        The reference sub's header, for :func:`write_stack` to build on.
    report : StackReport
    """
    import astroalign as aa
    from skimage.transform import warp

    paths = [Path(p) for p in paths]
    if not paths:
        raise ValueError("no sub-exposures given")

    acc = wacc = None
    ref_half = ref_header = None
    n_ok = n_bad = 0
    resids, rotations, reasons = [], [], []
    total_exptime = 0.0

    for k, path in enumerate(paths):
        with fits.open(path, memmap=False) as hdul:
            raw = hdul[0].data
            header = hdul[0].header
            exptime = header.get("EXPTIME", header.get("EXPOSURE")) or 0.0
            if k == 0:
                ref_header = header.copy()
                pattern = pattern or pattern_of(header)

        if k == 0:
            ny, nx = raw.shape
            acc = np.zeros((3, ny, nx), dtype=np.float64)
            wacc = np.zeros((ny, nx), dtype=np.float64)
            ref_half = green_half(raw, pattern)
            acc += debayer(raw, pattern)
            wacc += 1.0
            n_ok += 1
            total_exptime += float(exptime)
            continue

        try:
            tf, (src, dst) = aa.find_transform(
                green_half(raw, pattern), ref_half, detection_sigma=detection_sigma
            )
            resid = float(np.median(np.linalg.norm(tf(src) - dst, axis=1)))
            if len(src) < min_match or resid > max_resid_px:
                raise RuntimeError(f"poor fit: {len(src)} matches, {resid:.2f} px residual")
        except Exception as exc:
            n_bad += 1
            reasons.append(f"{path.name}: {type(exc).__name__}: {exc}")
            continue

        full = _scale_transform(tf)
        cube = debayer(raw, pattern)
        # warp maps output coordinates -> input coordinates, hence the inverse.
        for b in range(3):
            acc[b] += warp(cube[b], full.inverse, order=1, mode="constant",
                           cval=0.0, preserve_range=True)
        wacc += warp(np.ones(acc.shape[1:], dtype=np.float32), full.inverse, order=1,
                     mode="constant", cval=0.0, preserve_range=True)
        n_ok += 1
        total_exptime += float(exptime)
        resids.append(resid)
        rotations.append(np.degrees(tf.rotation))

    covered = wacc > MIN_WEIGHT
    cube = np.zeros_like(acc, dtype=np.float32)
    for b in range(3):
        cube[b] = np.where(covered, acc[b] / np.maximum(wacc, 1e-6), 0.0)

    report = StackReport(
        n_in=len(paths), n_ok=n_ok, n_bad=n_bad,
        resid_median_px=float(np.median(resids)) if resids else 0.0,
        rot_span_deg=float(np.ptp(rotations)) if rotations else 0.0,
        cover_frac=float(covered.mean()),
        total_exptime=total_exptime,
        reasons=reasons,
    )
    return cube, wacc, ref_header, report


def stack_frame(paths, **kwargs):
    """Co-add subs and return the result as a :class:`frames.SeestarFrame`.

    The convenient entry point: the result drops straight into
    :func:`photometry.extract_sources` and everything downstream, exactly like an on-board
    stack. Returns ``(frame, report)``.
    """
    from .frames import SeestarFrame, _model_of

    cube, _weight, header, report = coadd(paths, **kwargs)
    header = _stack_header(header, report)
    return SeestarFrame(
        data=cube, header=header, model=_model_of(header),
        path=Path(paths[0]), layout="stacked",
    ), report


def _stack_header(header, report):
    """The reference header, updated to describe the stack rather than one sub.

    ``STACKCNT`` and ``TOTALEXP`` are written in the *native Seestar dialect* so
    ``frames.frame_metadata`` reads a locally-built stack exactly as it reads an on-board
    one -- otherwise the per-sub ``EXPTIME`` would be mistaken for the total.
    """
    out = header.copy()
    out["STACKCNT"] = (report.n_ok, "subs successfully registered and co-added")
    out["TOTALEXP"] = (report.total_exptime, "total on-sky exposure [s]")
    out["STACKER"] = ("seestar-photometry", "local stack: astroalign + weighted mean")
    out["NSUBSIN"] = (report.n_in, "sub-exposures offered")
    out["NSUBSBAD"] = (report.n_bad, "subs rejected at registration")
    out["REGRESID"] = (round(report.resid_median_px, 4), "median registration residual [px]")
    out["ROTSPAN"] = (round(report.rot_span_deg, 4), "field rotation across the stack [deg]")
    out["COVERFRC"] = (round(report.cover_frac, 4), "fraction of pixels covered")
    # The reference sub's own WCS no longer describes the co-add's astrometry well enough
    # for photometry, and leaving it in place invites exactly the ~1 arcmin error the
    # package re-solves to avoid. Drop it so the frame reads as unsolved.
    for key in ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
                "CD1_1", "CD1_2", "CD2_1", "CD2_2", "CDELT1", "CDELT2", "PLTSOLVD"):
        out.pop(key, None)
    return out


def write_stack(path, cube, header, report=None, weight=None, overwrite=True):
    """Write a stack as a native-layout Seestar FITS cube.

    Output matches the on-board layout -- a single ``(3, ny, nx)`` primary HDU -- so it is
    indistinguishable to anything reading it. Pass ``weight`` to append the coverage map as
    a ``WEIGHT`` extension; the loader ignores it as a non-science plane.
    """
    header = _stack_header(header, report) if report is not None else header.copy()
    hdus = [fits.PrimaryHDU(data=np.asarray(cube, dtype=np.float32), header=header)]
    if weight is not None:
        hdus.append(fits.ImageHDU(data=np.asarray(weight, dtype=np.float32), name="WEIGHT"))
    fits.HDUList(hdus).writeto(path, overwrite=overwrite)
    return Path(path)
