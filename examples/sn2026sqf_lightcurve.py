#!/usr/bin/env python
"""SN 2026sqf in NGC 3310: a supernova on a bright host.

Differs from a field-star light curve in two ways, and both are handled by the package
rather than by bespoke code here:

* the frames are CrowdSky multi-extension files, already plate-solved server-side, so
  the WCS is **lifted** from the header rather than re-solved;
* the SN sits on the host galaxy, so host flux inside the aperture must be subtracted --
  empirically, at the SN's galactocentric radius, because NGC 3310's circumnuclear
  starburst ring sits almost exactly there and a smooth model under-predicts it.

There is no ensemble differential here: a supernova has no stable comparison of similar
brightness and the interest is the absolute magnitude, so each frame's own zero point is
used. Expected V ~ 12.95 near maximum.

Usage:
    uv run python examples/sn2026sqf_lightcurve.py [DATA_ROOT] [WORK_DIR]
"""

import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

from seestar_photometry import (
    LocalTree, Project, Target, astrometry, calibration, contamination, frames,
    photometry, pipeline, report,
)

# TNS / NED positions.
SN = Target("SN 2026sqf", ra=159.699838, dec=53.509472)   # 10h38m47.961s +53d30'34.10"
NUCLEUS = (159.6915, 53.503472)                            # NGC 3310 nucleus

DEFAULT_DATA = r"D:\seestar_paper_data_sets_2\crowdsky_ngc3310"
DEFAULT_WORK = r"D:\work\seestar_photometry\sn2026sqf"


def measure_epoch(proj, key, catalogue):
    """One frame -> one calibrated, host-corrected V magnitude."""
    from seestar_photometry import lightcurves

    frame = frames.load_frame(proj.source.path(key))
    wcs = astrometry.load_wcs(frame)
    if wcs is None:
        return None

    ext = photometry.extract_sources(
        frame, thresh=proj.thresh, enclosed=proj.enclosed_lightcurve
    )
    ext.match_gaia(catalogue, wcs=wcs, tol_arcsec=proj.match_tol_arcsec)
    cal = calibration.fit_zeropoint(ext.sources, band="G", mag_range=proj.fit_mag_range)

    aperture = float(ext.aperture[1])
    sn_xy = wcs.world_to_pixel_values(SN.ra, SN.dec)
    nuc_xy = wcs.world_to_pixel_values(*NUCLEUS)

    # Background-subtract the green plane once, then measure both the SN and the host.
    import sep

    green = np.asarray(frame.g, dtype=float)
    bkg = sep.Background(green)
    data_sub = green - bkg
    mask = contamination.star_mask(data_sub, bkg.globalrms, nuc_xy)
    host = contamination.galaxy_contamination(
        data_sub, nuc_xy, sn_xy, aperture, mask=mask
    )

    forced = photometry.forced_photometry(
        frame, [SN.ra], [SN.dec], wcs, enclosed=proj.enclosed_lightcurve
    )
    green_row = forced[np.asarray(forced["band"]) == "G"][0]
    flux = float(green_row["flux"])
    snr = float(green_row["snr"])

    corrected = flux - host["adu"]
    v_raw = float(photometry.instrumental_mag([flux])[0]) + cal.zeropoint
    v_corr = (cal.zeropoint - 2.5 * np.log10(corrected)) if corrected > 0 else np.nan

    meta = frames.frame_metadata(frame)
    times = lightcurves.frame_times(frame.header, _skycoord(), meta["total_exptime"])
    return {
        "frame": Path(str(frame.path)).name,
        "mjd_mid": times["mjd_mid"],
        "bjd_tdb": times["bjd_tdb"],
        "time_source": times["time_source"],
        "v_mag": v_corr,
        "v_raw": v_raw,
        # Photon noise on the SN, the zero-point error, and the azimuthal scatter of the
        # host estimate -- the last is a real systematic, not a rounding detail.
        "v_err": float(np.sqrt(
            (1.0857 / snr) ** 2
            + cal.zeropoint_err ** 2
            + (1.0857 * host["std"] / corrected) ** 2 if corrected > 0 else np.nan
        )),
        "v_err_sys": float(cal.rms),
        "host_adu": host["adu"],
        "host_frac": host["adu"] / flux if flux else np.nan,
        "n_azimuth": host["n_azimuth"],
        "zeropoint": cal.zeropoint,
        "n_cal": cal.n_stars,
        "total_exptime": meta["total_exptime"],
        "telescope": str(frame.header.get("TELESCOP", "")),
    }, (ext, cal, frame, data_sub, nuc_xy, sn_xy, aperture, host)


def _skycoord():
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    return SkyCoord(SN.ra * u.deg, SN.dec * u.deg)


def main():
    data_root = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_DATA
    work_dir = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_WORK

    proj = Project(
        target=SN,
        source=LocalTree(roots=[data_root]),
        work_dir=work_dir,
        # CrowdSky frames arrive plate-solved (PLTSOLVD), so lift rather than re-solve.
        solver="lift",
        catalogue_half_deg=1.0,
    )
    keys = proj.frames()
    if not keys:
        sys.exit(f"no frames under {data_root}")
    print(f"{len(keys)} frames | solver={proj.solver}")

    pipeline.solve_all(proj)
    catalogue = proj.catalogue()

    rows, latest = [], None
    for key in keys:
        try:
            result = measure_epoch(proj, key, catalogue)
        except Exception as exc:
            print(f"  failed {Path(str(key)).name}: {exc!r}")
            continue
        if result is None:
            continue
        row, artefacts = result
        rows.append(row)
        if latest is None or row["mjd_mid"] > latest[0]["mjd_mid"]:
            latest = (row, artefacts)
        print(f"  {row['frame'][:44]:46s} V={row['v_mag']:.3f} "
              f"host={row['host_frac']:.2%} ({row['n_azimuth']} azimuths)")

    if not rows:
        sys.exit("no epochs measured")

    table = Table(rows)
    table.sort("mjd_mid")
    out = proj.work_dir / "sn_lightcurve.ecsv"
    table.write(out, format="ascii.ecsv", overwrite=True)

    v = np.asarray(table["v_mag"], dtype=float)
    print(f"\n{len(table)} epochs | V {np.nanmin(v):.2f} to {np.nanmax(v):.2f} "
          f"(brightest {np.nanmin(v):.2f})")
    print("expected: V ~ 12.95 near maximum")
    print(f"host contamination: median {np.nanmedian(table['host_frac']):.1%} "
          "of the aperture flux")
    print(f"wrote {out}")

    # Diagnostics for the most recent epoch, plus the host cutout that shows exactly
    # what the contamination estimate sampled.
    if latest is not None:
        row, (ext, cal, frame, data_sub, nuc_xy, sn_xy, aperture, host) = latest
        report.frame_report(frame, ext, cal, proj.diagnostics_dir)
        cut = contamination.cutout(data_sub, nuc_xy, sn_xy, aperture)
        from seestar_photometry import plots

        ax = plots.host_cutout(cut, host)
        ax.figure.savefig(proj.diagnostics_dir / "sn_host_cutout.png", dpi=130,
                          bbox_inches="tight")
        print(f"figures in {proj.diagnostics_dir}")


if __name__ == "__main__":
    main()
