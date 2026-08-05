"""Project configuration: everything about a dataset that isn't code.

One :class:`Project` describes a target, where its frames live, where derived
products go, and the handful of numerical choices worth varying. Everything the old
scripts held in module-level globals lives here instead, so retargeting the pipeline
at a new object means constructing a different ``Project``, never editing code.

Both classes are thin data containers. The pipeline stages are module-level functions
in :mod:`pipeline` that take a project as their first argument.
"""

import os
from dataclasses import dataclass, field
from pathlib import Path

from .astrometry import ASTAP_EXE
from .calibration import FIT_MAG_RANGE


@dataclass
class Target:
    """The object a project is about.

    Parameters
    ----------
    name : str
        Human label, used in figure titles and filenames.
    ra, dec : float
        ICRS position in degrees. Used for the catalogue query centre, the
        barycentric time correction, and to identify the target in the stars table.
    source_id : int, optional
        Catalogue identifier, when known. Strongly preferred over letting the
        pipeline guess "the source nearest the pointing" -- that heuristic is right
        for a field centred on an isolated target and wrong in a crowded one.
    """

    name: str
    ra: float
    dec: float
    source_id: int = None

    @property
    def radec(self):
        return (self.ra, self.dec)


@dataclass
class Project:
    """A dataset to be reduced, and the choices that govern how.

    Parameters
    ----------
    target : Target
    source : FrameSource
        Where the frames come from -- e.g. :class:`frames.LocalTree`.
    work_dir : path-like
        Everything derived lands here: the catalogue cache, ``frames.ecsv``,
        ``stars.ecsv``, ``measurements.ecsv``, the error log and ``diagnostics/``.
        Never the data tree, and never the FITS headers. (The per-frame ``.wcs``
        sidecar is the one exception -- it lives beside its frame, because it is
        expensive to recompute and useful to every project that touches that frame.)
    solver : str
        ``"astap"`` (default, local and offline), ``"local"`` (anchored on the reference
        catalogue -- no binary, no network, no index files), ``"nova"``
        (astrometry.net), or ``"lift"`` for frames that already carry a trustworthy
        header solution.
    astap_exe : str
        Path to the ASTAP command-line binary.
    enclosed_characterise : float
        Enclosed-flux fraction sizing the aperture for frame characterisation.
        0.90 is the documented default and what the published depth numbers use.
    enclosed_lightcurve : float
        Enclosed-flux fraction for light-curve photometry. 0.95, because an aperture
        sweep showed the differential scatter floor minimises there: a slightly larger
        aperture damps the position-dependent aperture loss caused by Alt-Az field
        rotation, at negligible sky-noise cost for photon-dominated targets. Both
        values are deliberate -- see ``docs/photometry-design.md``.
    fit_mag_range : tuple
        Reference-V window for the zero-point fit: above saturation, below the
        detection floor.
    thresh : float
        SEP detection threshold in units of the background RMS. Kept low so the
        SNR > 5 cut, not detection, defines the sample.
    catalogue_half_deg : float
        Half-size of the cached catalogue box. Must cover the full dithered area, not
        just one FOV.
    catalogue_tiles : int
        Tiles per axis for the catalogue query. Leave at 1 unless the box is large;
        Gaia TAP is unreliable under concurrency.
    gmag_limit : float
        Faint limit of the catalogue query. Fainter costs query time and memory
        without adding usable calibrators.
    match_tol_arcsec : float
        Cross-match radius.
    catalogue_backend : str
        Where the reference catalogue comes from. ``"auto"`` (default) uses the offline
        copy when :mod:`gaiadb` has this field on disk and falls back to a Gaia TAP
        query when it does not, so installing the download is the only step needed to
        stop using the network. ``"local"`` and ``"tap"`` force one or the other;
        ``"local"`` raises rather than querying, which is what you want on a machine
        that must not reach the internet.
    epoch : float
        Decimal year to propagate catalogue positions to, e.g. ``2026.4``. Gaia
        positions are J2016.0, and a decade of proper motion moves the fastest stars by
        more than ``match_tol_arcsec``. Only the local backend can do this -- the TAP
        query does not fetch proper motions -- so it is ignored with a TAP catalogue.
        ``None`` leaves positions at the Gaia epoch.
    provenance : callable, optional
        ``provenance(frame) -> dict`` of extra frame-table columns (dataset name,
        stacking manifest fields, ...). Keeps project bookkeeping out of the core
        schema.
    """

    target: Target
    source: object
    work_dir: str
    solver: str = "astap"
    astap_exe: str = ASTAP_EXE
    enclosed_characterise: float = 0.90
    enclosed_lightcurve: float = 0.95
    fit_mag_range: tuple = FIT_MAG_RANGE
    thresh: float = 2.0
    catalogue_half_deg: float = 1.5
    catalogue_tiles: int = 1
    gmag_limit: float = 17.0
    match_tol_arcsec: float = 2.0
    catalogue_backend: str = "auto"
    epoch: float = None
    provenance: object = field(default=None, repr=False)

    def __post_init__(self):
        self.work_dir = Path(self.work_dir)

    # --- derived paths ---------------------------------------------------------------

    @property
    def catalogue_path(self):
        """ECSV cache of the reference catalogue for this field."""
        return self.work_dir / f"catalogue_{_slug(self.target.name)}.ecsv"

    @property
    def frames_path(self):
        """Per-frame quality table."""
        return self.work_dir / "frames.ecsv"

    @property
    def stars_path(self):
        """One row per measured catalogue source."""
        return self.work_dir / "stars.ecsv"

    @property
    def measurements_path(self):
        """Long table, one row per (source, frame, band)."""
        return self.work_dir / "measurements.ecsv"

    @property
    def diagnostics_dir(self):
        """Where diagnostic figures are saved."""
        return self.work_dir / "diagnostics"

    @property
    def log_path(self):
        return self.work_dir / "errors.log"

    def ensure_dirs(self):
        """Create ``work_dir`` (and ``diagnostics/``) if they don't exist."""
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.diagnostics_dir.mkdir(parents=True, exist_ok=True)
        return self

    # --- convenience ------------------------------------------------------------------

    def frames(self):
        """Keys of every frame in the source, in deterministic order."""
        return list(self.source.keys())

    @property
    def catalogue_radius_deg(self):
        """Cone radius covering the whole catalogue box, corners included.

        The same figure :func:`catalogs.fetch_gaia_mosaic` computes for a single tile,
        so the two backends cover the same sky and can be compared row for row.
        """
        return self.catalogue_half_deg * 2 ** 0.5 + 0.1

    def catalogue_backend_used(self):
        """Which backend :meth:`catalogue` would use. Touches no network."""
        if self.catalogue_backend != "auto":
            return self.catalogue_backend
        from . import gaiadb

        return "local" if gaiadb.covers(
            self.target.radec, self.catalogue_radius_deg) else "tap"

    def catalogue(self, overwrite=False):
        """Load (or build once) the cached reference catalogue for this field.

        Either backend writes the same ECSV to :attr:`catalogue_path`, so everything
        downstream is unaffected by which one ran.
        """
        from . import catalogs

        self.ensure_dirs()
        if self.catalogue_path.exists() and not overwrite:
            return catalogs.load_catalogue(self.catalogue_path)

        if self.catalogue_backend_used() == "local":
            from . import gaiadb

            catalogue = gaiadb.cone(
                self.target.radec, self.catalogue_radius_deg,
                gmag_limit=self.gmag_limit, epoch=self.epoch,
            )
            catalogue.write(self.catalogue_path, format="ascii.ecsv", overwrite=True)
            return catalogue

        return catalogs.fetch_gaia_mosaic(
            self.target.radec, self.catalogue_path,
            half_size_deg=self.catalogue_half_deg, n_tiles=self.catalogue_tiles,
            gmag_limit=self.gmag_limit, overwrite=overwrite,
        )

    @property
    def api_key(self):
        """astrometry.net key from the environment (only needed for ``solver="nova"``)."""
        return os.environ.get("ASTROMETRY_KEY")


def _slug(text):
    """Filename-safe lowercase slug, e.g. ``"MW Cam"`` -> ``"mw_cam"``."""
    keep = [c if (c.isalnum() or c in "-_") else "_" for c in str(text).strip().lower()]
    return "".join(keep).strip("_") or "field"
