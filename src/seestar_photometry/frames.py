"""Loading Seestar frames, and discovering them on disk.

Two FITS layouts reach this package, and both normalise to one in-memory
representation -- a ``(3, ny, nx)`` float32 array in R, G, B order:

``"cube"`` -- **native Seestar**
    The telescope writes the debayered, stacked RGB planes as a single 3-D image in
    the primary HDU: FITS axes ``(NAXIS1=nx, NAXIS2=ny, NAXIS3=3)``, which astropy
    reads into numpy as ``(3, ny, nx)``. Stored uint16. Frame size depends on the
    model and binning (S50: 1080x1920; S30pro: 2160x3840, or 1080x1920 binned).

``"mef"`` -- **CrowdSky**
    The CrowdSky platform re-stacks and plate-solves, then writes a multi-extension
    file: an **empty** primary HDU carrying all the metadata, then ``RED``,
    ``GREEN``, ``BLUE`` ``ImageHDU``\\ s, usually a ``FOOTPRINT`` coverage plane, and
    a ``STAR-TAB`` binary table holding the server's own SEP catalogue already
    cross-matched to Gaia.

Because the ``"mef"`` primary HDU holds no image, **the image shape must always come
from the data, never from a header**.

The two layouts also speak different header dialects for exposure metadata, and
``EXPTIME`` means different things in each -- see :func:`frame_metadata`. Read
metadata through that function, never straight off the header.
"""

import glob
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol, runtime_checkable

import numpy as np
from astropy.io import fits

#: Order of the colour planes along axis 0. Green (index 1) is the science band.
BANDS = ("R", "G", "B")

#: ``EXTNAME``\\ s accepted for each science plane in the multi-extension layout,
#: in ``BANDS`` order. Matched case-insensitively.
_PLANE_NAMES = (("RED", "R"), ("GREEN", "G"), ("BLUE", "B"))

#: Extension holding the server's own source catalogue (multi-extension layout).
_STAR_TAB = "STAR-TAB"

#: Derived planes that must never be mistaken for a science plane.
_NON_SCIENCE = ("FOOTPRINT", "MASK", "WEIGHT")


@dataclass
class SeestarFrame:
    """A single stacked Seestar exposure.

    A thin data container. The analysis lives in :mod:`photometry`,
    :mod:`astrometry` and friends, which all take a frame as their first argument.

    Attributes
    ----------
    data : np.ndarray
        The ``(3, ny, nx)`` R, G, B cube as float32, whatever the file layout.
    header : fits.Header
        The primary HDU header -- the metadata carrier in both layouts.
    model : str
        Seestar model, ``"S50"`` or ``"S30pro"``.
    path : Path
        Path the frame was loaded from.
    layout : str
        ``"cube"`` (native) or ``"mef"`` (CrowdSky). See the module docstring.
    star_tab : astropy.table.Table or None
        The ``STAR-TAB`` catalogue, when the file carries one (CrowdSky only): the
        server's SEP detections with ``x``, ``y``, ``flux``, ``ra``, ``dec`` and
        Gaia cross-match columns. Useful as an independent cross-check of our own
        extraction, and as a source list for a WCS solve.
    """

    data: np.ndarray
    header: fits.Header
    model: str
    path: Path
    layout: str = "cube"
    star_tab: object = field(default=None, repr=False)

    @property
    def r(self):
        return self.data[0]

    @property
    def g(self):
        """The green plane -- the science band."""
        return self.data[1]

    @property
    def b(self):
        return self.data[2]

    @property
    def shape(self):
        """Image shape ``(ny, nx)``, from the data (headers can't be trusted)."""
        return self.data.shape[1:]


def _normalise_cube(data):
    """Return a ``(3, ny, nx)`` float32 cube from a 3-D array of either ordering.

    Every Seestar file observed so far is channel-*first* once astropy has read it
    (FITS ``NAXIS3=3`` becomes numpy axis 0), and that is the fast path. The
    channel-*last* transpose is cheap insurance for exports that have been through
    an image library, where ``(ny, nx, 3)`` is the natural convention -- getting
    this wrong silently swaps a colour plane for an image row, which would look
    like a catastrophic PSF rather than a load error.
    """
    if data.ndim != 3:
        raise ValueError(f"expected a 3-D RGB cube, got shape {data.shape}")
    if data.shape[0] == 3:
        cube = data
    elif data.shape[-1] == 3:
        cube = np.moveaxis(data, -1, 0)
    else:
        raise ValueError(f"no length-3 colour axis in cube of shape {data.shape}")
    return np.ascontiguousarray(cube, dtype=np.float32)


def _plane_index(hdul):
    """Map each band to its HDU index in a multi-extension file.

    Resolved by ``EXTNAME`` (``RED``/``GREEN``/``BLUE``, or bare ``R``/``G``/``B``),
    falling back to the first three 2-D image HDUs in file order. ``FOOTPRINT`` and
    other derived planes are never candidates -- including one would shift the
    colour assignment by one plane.
    """
    by_name = {}
    for i, hdu in enumerate(hdul):
        name = str(hdu.header.get("EXTNAME", "")).upper()
        if name and name not in _NON_SCIENCE:
            by_name.setdefault(name, i)

    found = []
    for aliases in _PLANE_NAMES:
        hit = next((by_name[a] for a in aliases if a in by_name), None)
        if hit is None:
            break
        found.append(hit)
    if len(found) == 3:
        return found

    candidates = [
        i for i, hdu in enumerate(hdul)
        if getattr(hdu, "data", None) is not None
        and len(getattr(hdu, "shape", ()) or ()) == 2
        and str(hdu.header.get("EXTNAME", "")).upper() not in _NON_SCIENCE
    ]
    if len(candidates) < 3:
        raise ValueError(
            f"need 3 science planes, found {len(candidates)} 2-D image HDUs"
        )
    return candidates[:3]


def _read_star_tab(hdul):
    """The ``STAR-TAB`` catalogue as a Table, or ``None`` if absent."""
    for hdu in hdul:
        if str(hdu.header.get("EXTNAME", "")).upper() == _STAR_TAB:
            from astropy.table import Table

            return Table(hdu.data)
    return None


def load_frame(path):
    """Load a Seestar frame from either FITS layout.

    The colour planes are promoted to float32 so downstream arithmetic can't
    overflow the native uint16 storage, and are always returned as
    ``(3, ny, nx)`` in R, G, B order. The header is always the primary HDU's.
    """
    path = Path(path)
    with fits.open(path) as hdul:
        header = hdul[0].header
        primary = hdul[0].data
        if primary is not None and primary.ndim == 3:
            layout = "cube"
            data = _normalise_cube(primary)
            star_tab = _read_star_tab(hdul)
        else:
            layout = "mef"
            idx = _plane_index(hdul)
            data = np.ascontiguousarray(
                np.stack([hdul[i].data for i in idx]), dtype=np.float32
            )
            star_tab = _read_star_tab(hdul)
    return SeestarFrame(
        data=data, header=header, model=_model_of(header), path=path,
        layout=layout, star_tab=star_tab,
    )


def _model_of(header):
    """Seestar model from the header.

    Keyed off ``TELESCOP`` (e.g. ``"S50_8a95aa90"``, the model plus the unit's
    serial), not ``INSTRUME``, because ``INSTRUME`` is inconsistent across firmware
    versions and absent from some exports.
    """
    telescop = str(header.get("TELESCOP", ""))
    return "S50" if "S50" in telescop else "S30pro"


def unit_id(header):
    """The individual telescope's serial, e.g. ``"8a95aa90"`` from ``TELESCOP``.

    Each physical Seestar has its own zero point, PSF and condition response, so
    per-unit grouping matters when combining data from several telescopes.
    Returns ``""`` if the header carries no serial.
    """
    telescop = str(header.get("TELESCOP", ""))
    return telescop.split("_", 1)[1] if "_" in telescop else ""


# --- exposure metadata: two header dialects ------------------------------------------

def _exposure(h):
    """Resolve ``(n_exp, exptime_per_sub, total_exptime)`` across header dialects.

    The two writers disagree, and ``EXPTIME`` is the trap -- it is the **per-sub**
    exposure in a native file and the **total on-sky** time in a CrowdSky one:

    ==================  ====================  ==========================
    quantity            native Seestar        CrowdSky
    ==================  ====================  ==========================
    sub-exposure count  ``STACKCNT``          ``NIMAGES``
    per-sub exposure    ``EXPTIME``           ``EXPTIME / NIMAGES``
    total on-sky        ``TOTALEXP``          ``EXPTIME``
    ==================  ====================  ==========================

    Taking ``EXPTIME`` at face value on a CrowdSky frame gives a 410 s sub-exposure
    and (via ``n_exp * exptime``) a ~4.6 hour "total" for a 15-minute stack. The
    dialect is decided by which keywords are present, not by the file layout, so a
    future exporter that mixes them still resolves correctly.
    """
    exptime_kw = h.get("EXPTIME", h.get("EXPOSURE"))
    if h.get("TOTALEXP") is not None:                       # native dialect
        n_exp, total, per_sub = h.get("STACKCNT"), h.get("TOTALEXP"), exptime_kw
    elif h.get("NIMAGES") is not None:                      # CrowdSky dialect
        n_exp, total = h.get("NIMAGES"), exptime_kw
        per_sub = (total / n_exp) if (total and n_exp) else None
    else:                                                   # sparse header
        n_exp, per_sub = h.get("STACKCNT", h.get("NIMAGES")), exptime_kw
        total = n_exp * per_sub if (n_exp is not None and per_sub is not None) else None
    return n_exp, per_sub, total


def exposure_span(header):
    """``(start, end)`` of the on-sky exposure as ISO strings, ``end`` may be None.

    CrowdSky records the true wall-clock span of the co-add in ``OB-START`` /
    ``OB-END``, which is strictly better than inferring a mid-exposure time from a
    single epoch plus the on-sky integration -- it includes the inter-sub overhead
    that is otherwise unknowable. Native frames carry only ``DATE-OBS``.
    """
    start = header.get("OB-START") or header.get("DATE-OBS")
    end = header.get("OB-END")
    norm = lambda v: str(v).replace(" ", "T") if v is not None else None
    return norm(start), norm(end)


def _airmass(header):
    """Airmass (sec z) from the header pointing, site and time; nan if unavailable.

    Computed rather than read -- Seestar headers carry no ``AIRMASS``. IERS
    auto-download is disabled and a stale bundled table is extrapolated without
    complaint; that is far more accuracy than sec z needs.
    """
    lat, lon = header.get("SITELAT"), header.get("SITELONG")
    if lat is None or lon is None or "DATE-OBS" not in header:
        return float("nan")
    import warnings

    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord
    from astropy.time import Time
    from astropy.utils import iers

    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        altaz = SkyCoord(header["RA"] * u.deg, header["DEC"] * u.deg).transform_to(
            AltAz(
                obstime=Time(str(header["DATE-OBS"]).replace(" ", "T")),
                location=EarthLocation(lat=lat * u.deg, lon=lon * u.deg),
            )
        )
        return float(altaz.secz)


def frame_metadata(frame):
    """Observing metadata from a frame's header, dialect-resolved.

    Returns a dict with ``n_exp`` (stacked sub-exposures), ``exptime`` (per sub, s),
    ``total_exptime`` (on-sky, s), ``obs_start``/``obs_end`` (ISO, ``obs_end`` may be
    None), ``airmass`` (computed), ``site_lat``/``site_lon`` (deg), ``ccd_temp``,
    ``eqmode`` (0 = Alt-Az, 1 = equatorial; absent from CrowdSky headers) and
    ``unit`` (the telescope serial).

    Missing values come back as ``nan`` (or ``None`` for the ISO strings) so sparser
    headers from other users' units still flow through.
    """
    h = frame.header
    n_exp, per_sub, total = _exposure(h)
    start, end = exposure_span(h)
    meta = {
        "n_exp": n_exp,
        "exptime": per_sub,
        "total_exptime": total,
        "airmass": _airmass(h),
        "site_lat": h.get("SITELAT"),
        "site_lon": h.get("SITELONG"),
        "ccd_temp": h.get("CCD-TEMP"),
        "eqmode": h.get("EQMODE"),  # 0 = Alt-Az, 1 = equatorial
    }
    meta = {k: (np.nan if v is None else v) for k, v in meta.items()}
    meta["obs_start"] = start
    meta["obs_end"] = end
    meta["unit"] = unit_id(h)
    return meta


# --- frame discovery -----------------------------------------------------------------

@runtime_checkable
class FrameSource(Protocol):
    """Where a project's frames come from.

    The pipeline only needs to enumerate frames and get a local path for each, so
    that is the whole protocol. It exists so a remote archive (e.g. CrowdSky, which
    discovers frames through an API and downloads them on demand) can be dropped in
    later without touching any pipeline code.

    Implementations must be picklable -- the batch runner ships them to worker
    processes.
    """

    def keys(self):
        """Stable identifiers for every available frame, in a deterministic order.

        A key must be usable as a filename component and must not change between
        runs; it is what makes the pipeline resumable.
        """

    def path(self, key):
        """A local filesystem path for one frame, materialising it if needed."""


@dataclass
class LocalTree:
    """Frames already on disk, discovered by walking one or more directory trees.

    Parameters
    ----------
    roots : sequence of path-like
        Directories to search recursively.
    patterns : tuple of str
        Filename globs to match. Defaults to Seestar's ``.fit`` and the more
        conventional ``.fits``.
    curate : callable, optional
        ``curate(path) -> bool``; frames for which it returns False are dropped.
        Datasets accumulate misfiled frames (a different target left in a folder)
        and deliberately-excluded ones, and that judgement is dataset-specific, so
        it is a hook rather than a built-in rule. Called once per frame during
        discovery, so keep it cheap -- read a header, not the pixels.
    """

    roots: tuple
    patterns: tuple = ("*.fit", "*.fits")
    curate: object = None

    def __post_init__(self):
        if isinstance(self.roots, (str, os.PathLike)):
            self.roots = (self.roots,)
        self.roots = tuple(self.roots)

    def keys(self):
        """Absolute paths of every matching frame, sorted, de-duplicated.

        The path *is* the key for a local tree. Sorting is case-insensitive so the
        order is stable across platforms; de-duplication guards against overlapping
        or nested roots quietly measuring a frame twice.
        """
        seen, out = set(), []
        for root in self.roots:
            for pattern in self.patterns:
                hits = glob.glob(os.path.join(str(root), "**", pattern), recursive=True)
                for p in hits:
                    key = os.path.normcase(os.path.abspath(p))
                    if key not in seen:
                        seen.add(key)
                        out.append(os.path.abspath(p))
        out.sort(key=str.lower)
        if self.curate is not None:
            out = [p for p in out if self.curate(p)]
        return out

    def path(self, key):
        """A local tree's frames are already local, so the key is the path."""
        return key


def object_name(path):
    """Normalised FITS ``OBJECT`` (upper case, spaces stripped), ``""`` if unreadable.

    A helper for writing :class:`LocalTree` curation predicates.
    """
    try:
        return str(fits.getheader(str(path)).get("OBJECT", "")).upper().replace(" ", "")
    except Exception:
        return ""
