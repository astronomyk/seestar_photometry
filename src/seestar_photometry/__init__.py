"""Time-domain photometry from ZWO Seestar smart-telescope stacks.

Calibrated onto Gaia DR3 synthetic Johnson V, with ensemble differential photometry
good to well under 1% on a well-chosen comparison ensemble.

A minimal session::

    from seestar_photometry import Project, Target, LocalTree, pipeline, lightcurves

    proj = Project(
        target   = Target("MW Cam", ra=186.6821, dec=81.474),
        source   = LocalTree(roots=[r"D:\\data\\MW Cam s50\\stacks"]),
        work_dir = r"D:\\work\\mwcam",
    )
    pipeline.solve_all(proj)                                   # .wcs sidecars
    frames = pipeline.build_frame_table(proj, diagnostics=3)    # frames.ecsv + figures
    stars, meas = pipeline.build_measurements(proj)             # the two LC tables

    comps = lightcurves.select_comparisons(stars, dmag=1.0, colour_tol=0.3)
    lc = lightcurves.differential_lightcurve(meas, lightcurves.target_id_of(stars), comps)

See ``CLAUDE.md`` for the conventions and ``docs/`` for why each numerical choice is
what it is.
"""

from . import (
    astrometry,
    calibration,
    catalogs,
    contamination,
    debayer,
    depth,
    examples,
    frames,
    gaiadb,
    kron,
    lightcurves,
    photometry,
    pipeline,
    plots,
    project,
    quality,
    report,
    stacking,
)
from .frames import BANDS, LocalTree, SeestarFrame, load_frame
from .project import Project, Target

try:
    from importlib.metadata import version as _pkg_version

    #: Read from installed metadata, so pyproject.toml is the single source of truth.
    #: A hardcoded literal here silently drifted out of sync once already.
    __version__ = _pkg_version("seestar-photometry")
except Exception:  # running from a source tree with no install
    __version__ = "0.0.0.dev0"

__all__ = [
    # submodules
    "astrometry", "calibration", "catalogs", "contamination", "debayer", "depth",
    "examples", "frames", "gaiadb", "kron", "lightcurves", "photometry", "pipeline",
    "plots", "project", "quality", "report", "stacking",
    # the names you actually reach for
    "BANDS", "LocalTree", "Project", "SeestarFrame", "Target", "load_frame",
]
