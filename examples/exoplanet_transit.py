#!/usr/bin/env python
"""Template: a new target from scratch. Edit the CONFIG block, run, read the figures.

Written around a transit, because that is the case with the tightest requirements --
sub-1% precision on a single night, with correct absolute timing. Nothing here is
transit-specific except the final plot and the timing emphasis; the same driver works
for any variable.

Usage:
    uv run python examples/exoplanet_transit.py
"""

import sys

import numpy as np

from seestar_photometry import (
    LocalTree, Project, Target, lightcurves, pipeline, report,
)

# ---------------------------------------------------------------- CONFIG: edit this
TARGET = Target(
    name="HD 189733",
    ra=300.18216,      # degrees, ICRS
    dec=22.71100,
    # source_id=1234567890123456789,   # set this if you know it -- see below
)

DATA_ROOTS = [r"D:\seestar_data\hd189733\stacks"]
WORK_DIR = r"D:\work\seestar_photometry\hd189733"

# The catalogue box must cover the whole dithered area, not one FOV. 1.5 deg half-size
# comfortably covers a night's dithering; enlarge only if the pointing wandered.
CATALOGUE_HALF_DEG = 1.5

# Comparison-ensemble cuts. Start here and tighten based on lc_comparison_grid.
COMP_DMAG = 1.0          # within +-1 mag of the target
COMP_COLOUR_TOL = 0.3    # within +-0.3 in BP-RP
COMP_MAX_SEP_ARCMIN = 15
# -----------------------------------------------------------------------------------


def main():
    proj = Project(
        target=TARGET,
        source=LocalTree(roots=DATA_ROOTS),
        work_dir=WORK_DIR,
        catalogue_half_deg=CATALOGUE_HALF_DEG,
        solver="astap",
    )

    keys = proj.frames()
    if not keys:
        sys.exit(f"no .fit/.fits frames found under {DATA_ROOTS}")
    print(f"{len(keys)} frames found")

    # Stage 1: plate solve. The on-board WCS is off by ~1 arcmin and unusable, so this
    # is not optional. Idempotent -- safe to re-run, safe to interrupt.
    pipeline.solve_all(proj)

    # Stage 2: characterise every frame. Read frames_summary.png before going further:
    # it tells you how much of the night was photometric.
    frames = pipeline.build_frame_table(proj, diagnostics=3)
    usable = int((np.asarray(frames["rms"], dtype=float) < 0.06).sum())
    print(f"\n{usable}/{len(frames)} frames are photometric-grade (rms < 0.06)")
    if usable < len(frames) // 2:
        print("  more than half the night is marginal -- check frames_zp_vs_time.png "
              "for cloud before trusting the light curve")

    # Stage 3: forced photometry of every catalogue source.
    stars, meas = pipeline.build_measurements(proj)

    # If you did not set Target.source_id, check the flagged target is really yours:
    # 'nearest to the pointing' is a heuristic, not a guarantee.
    target = stars[np.asarray(stars["is_target"])]
    print(f"\ntarget: source_id={int(target['source_id'][0])} "
          f"V={float(target['v_jkc_mag'][0]):.2f} "
          f"at {float(target['sep_target_arcmin'][0]):.2f}' from the pointing")

    comps = lightcurves.select_comparisons(
        stars, dmag=COMP_DMAG, colour_tol=COMP_COLOUR_TOL,
        max_sep_arcmin=COMP_MAX_SEP_ARCMIN,
    )
    print(f"{len(comps)} comparison stars")
    if len(comps) < 5:
        print("  few comparisons -- loosen COMP_DMAG or COMP_MAX_SEP_ARCMIN")

    lc = lightcurves.differential_lightcurve(
        meas, lightcurves.target_id_of(stars), comps, band="G",
    )
    if not len(lc):
        sys.exit("no epochs survived; check errors.log and the diagnostics")

    scatter = lc.meta["scatter"] * 1000
    median_err = float(np.nanmedian(np.asarray(lc["mag_err"], dtype=float))) * 1000
    print(f"\n{len(lc)} epochs | scatter {scatter:.1f} mmag | "
          f"median error {median_err:.1f} mmag")
    if scatter > 2 * median_err:
        print("  scatter exceeds the photon prediction -- a systematic remains. "
              "Read lc_comparison_grid.png: a bad comparison is the usual cause.")

    # For a transit the timing is the point, so state which timing path was taken.
    span_hours = (float(np.max(lc["time"])) - float(np.min(lc["time"]))) * 24.0
    print(f"baseline {span_hours:.2f} h, BJD_TDB "
          f"{float(np.min(lc['time'])):.5f} to {float(np.max(lc['time'])):.5f}")

    report.lightcurve_report(
        lc, stars, meas, comps, proj.diagnostics_dir, project=proj,
    )
    lc.write(proj.work_dir / "lightcurve.ecsv", format="ascii.ecsv", overwrite=True)
    print(f"\nwrote {proj.work_dir / 'lightcurve.ecsv'}")
    print(f"figures in {proj.diagnostics_dir} -- read lc_comparison_grid.png first")


if __name__ == "__main__":
    main()
