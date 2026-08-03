#!/usr/bin/env python
"""MW Cam: a delta-Scuti light curve from Seestar S50 stacks.

The reference driver, and the parity check against the predecessor package. The expected
result on the c17 Alt-Az set is P = 0.1294 d with a scatter floor of ~23 mmag.

Everything dataset-specific is in ``Project`` below -- nothing else in this file needs
editing to point it at a different target.

Usage:
    uv run python examples/mwcam_lightcurve.py [DATA_ROOT] [WORK_DIR]
"""

import sys
from pathlib import Path

from seestar_photometry import (
    LocalTree, Project, Target, lightcurves, pipeline, report,
)

# MW Cam (RA, Dec) in degrees, and its Gaia DR3 identifier. Naming the source_id is
# better than letting the pipeline pick "nearest to the pointing" -- that heuristic is
# right for an isolated centred target and wrong in a crowded field.
TARGET = Target("MW Cam", ra=186.6821, dec=81.474)

DEFAULT_DATA = r"D:\seestar_paper_data_sets_2\MW Cam s50\stacks"
DEFAULT_WORK = r"D:\work\seestar_photometry\mwcam"


def only_mw_cam(path):
    """Drop frames whose OBJECT positively names a different field.

    A blocklist, not an allowlist, and deliberately so: MW Cam labels are inconsistent
    across contributors ('MW Cam', 'MWCam_flat_NN', 'mw', 'mw correct') and some units
    write no OBJECT at all. Anything MW-ish or blank is kept; only frames that clearly
    name another target are dropped -- in practice a misfiled IC 3568 set and a stray
    NGC 3310 frame.
    """
    from seestar_photometry.frames import object_name

    obj = object_name(path)
    return obj == "" or "MW" in obj


def main():
    data_root = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_DATA
    work_dir = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_WORK

    proj = Project(
        target=TARGET,
        source=LocalTree(roots=[data_root], curate=only_mw_cam),
        work_dir=work_dir,
        solver="astap",
    )
    print(f"{len(proj.frames())} frames | work_dir {proj.work_dir}")

    # 1. WCS sidecars. Idempotent -- already-solved frames are skipped.
    pipeline.solve_all(proj)

    # 2. Per-frame characterisation, with diagnostics for the first few frames.
    frames = pipeline.build_frame_table(proj, diagnostics=3)

    # 3. Forced photometry of every catalogue source in every frame.
    stars, meas = pipeline.build_measurements(proj)

    # 4. A comparison ensemble: close, similar brightness, similar colour, not variable.
    comps = lightcurves.select_comparisons(
        stars, dmag=1.0, colour_tol=0.3, max_sep_arcmin=15,
    )
    target_id = lightcurves.target_id_of(stars)
    print(f"target source_id={target_id} | {len(comps)} comparisons")

    lc = lightcurves.differential_lightcurve(meas, target_id, comps, band="G")
    if not len(lc):
        sys.exit("no epochs survived; check errors.log and the diagnostics")

    pg = lightcurves.periodogram(lc, min_period=0.02, max_period=1.0)
    print(f"\n{len(lc)} epochs | scatter {lc.meta['scatter'] * 1000:.1f} mmag")
    print(f"best period {pg['best_period']:.5f} d | FAP {pg['fap']:.2e}")
    print("expected: P = 0.1294 d, scatter ~23 mmag")

    # 5. The figures. Read lc_comparison_grid first -- it is the check that matters.
    reference = _reference_frame(proj)
    report.lightcurve_report(
        lc, stars, meas, comps, proj.diagnostics_dir, project=proj,
        frame=reference[0], wcs=reference[1],
    )
    lc.write(proj.work_dir / "lightcurve.ecsv", format="ascii.ecsv", overwrite=True)
    print(f"\nwrote {proj.work_dir / 'lightcurve.ecsv'}")
    print(f"figures in {proj.diagnostics_dir}")


def _reference_frame(proj):
    """The first solved frame and its WCS, for the finder chart."""
    from seestar_photometry import astrometry, frames as frames_mod

    for key in proj.frames():
        frame = frames_mod.load_frame(proj.source.path(key))
        wcs = astrometry.load_wcs(frame)
        if wcs is not None:
            return frame, wcs
    return None, None


if __name__ == "__main__":
    main()
