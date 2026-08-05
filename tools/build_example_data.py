#!/usr/bin/env python
"""Regenerate the example dataset in ``example_data/`` and pack it for release.

Maintenance script, not shipped in the wheel. It needs the full MW Cam datasets, which
live outside this repo; run it only when the example data must change.

What it produces, and why each piece is there:

``stack_c17_15min``   a clean 760 s c17 stack (rms 0.016 mag, 59 calibration stars)
                      -- the workhorse for photometry, zero point and light-curve docs
``stack_c17_30min``   a 1460 s stack of the same field, two hours later
                      -- gives depth-vs-exposure and a two-epoch light curve something
                         real to compare against
``stack_saturated``   a 280 s Zcom20 stack whose brightest star clips
                      -- the saturation-limit use case needs a frame that saturates
``crowdsky_mef``      one CrowdSky multi-extension frame, server plate-solved
                      -- exercises the second FITS layout and the lifted-WCS path
``raw_sub_1..5``      five consecutive 20 s c17 raw subs, same night as the stacks
                      -- the debayer and local-stacking use case
``gaia_mwcam``        Gaia DR3 rows covering the cutout footprint, with synthetic JKC V
                      -- makes zero-point calibration work offline, with no TAP query
``*.wcs``             per-frame sidecars, CRPIX shifted to the cutout

Everything is a **1000x1000 cutout** of the real frame, gzipped. Real pixels, real
headers, real stars; the WCS stays valid because the cutout origin is subtracted from
CRPIX. Full frames would be ~15 MB each.

The dataset is **not** shipped in the wheel and **not** committed: this script writes it to
``example_data/`` and packs it into ``example-data-<version>.tar.gz``, which is published as
a GitHub release asset and fetched on demand by ``examples.download()``. After running this,
paste the printed SHA-256 into ``examples.DATA_SHA256`` and upload the archive:

    gh release create example-data-v1 example_data/example-data-v1.tar.gz         --title "Example data v1" --notes "Real Seestar cutouts for tests and docs."

Usage:
    uv run python tools/build_example_data.py [--size 1000]
"""

from __future__ import annotations

import argparse
import glob
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "example_data"

DATA_ROOT = Path(r"D:\seestar_paper_data_sets_2")
S50 = DATA_ROOT / "MW Cam s50" / "stacks"
RAWS = DATA_ROOT / "MW Cam s50" / "raws" / "MW Cam_sub s50 c17"
CROWDSKY = DATA_ROOT / "crowdsky_mwcam"
GAIA = Path(r"D:\tmp\photometry_adventures_2\paper\data\gaia_MWCam.ecsv")

#: Chosen from the per-frame quality table -- see the module docstring.
STACKS = {
    "stack_c17_15min": "c17/**/c17_X015min_20260503-0045_n38_t760s_20s_IRCUT.fit",
    "stack_c17_30min": "c17/**/c17_X030min_20260503-0130_n73_t1460s_20s_IRCUT.fit",
    "stack_saturated": "Zcom20/**/Zcom20_X005min_20260627-0050_n14_t280s_20s_IRCUT.fit",
}
N_RAWS = 5

#: Must match ``examples.DATA_VERSION``.
DATA_VERSION = "v1"

#: Gaia columns worth shipping. Dropping the rest keeps the table small.
GAIA_COLUMNS = ("source_id", "ra", "dec", "phot_g_mean_mag", "phot_bp_mean_mag",
                "phot_rp_mean_mag", "bp_rp", "phot_variable_flag",
                "v_jkc_mag", "b_jkc_mag", "r_jkc_mag")


def centre_slice(shape, size, centre=None):
    """``(y0, y1, x0, x1)`` for a ``size x size`` cutout, clipped inside the frame.

    Centred on the frame unless ``centre=(y, x)`` is given.
    """
    ny, nx = shape
    if centre is None:
        y0, x0 = (ny - size) // 2, (nx - size) // 2
    else:
        y0, x0 = int(centre[0]) - size // 2, int(centre[1]) - size // 2
    y0 = int(np.clip(y0, 0, max(ny - size, 0)))
    x0 = int(np.clip(x0, 0, max(nx - size, 0)))
    return y0, min(y0 + size, ny), x0, min(x0 + size, nx)


def brightest_pixel(cube):
    """``(y, x)`` of the peak of the green plane.

    Used to place the saturation-example cutout. Centring on the frame is right for every
    other frame -- it keeps the same patch of sky across them -- but the star that clips
    the sensor is wherever it happens to be, and a cutout that misses it makes
    ``saturation_mag`` return nan, which is exactly the thing the example must show.
    """
    green = cube[1] if cube.ndim == 3 else cube
    return np.unravel_index(int(np.argmax(green)), green.shape)


def shift_wcs(header, x0, y0):
    """Move a header's WCS reference pixel to the cutout's coordinate origin.

    FITS CRPIX is 1-indexed, and the cutout starts at 0-indexed ``(x0, y0)``, so the new
    reference pixel is ``CRPIX - (x0, y0)``. Nothing else in the WCS changes: the plate
    scale, rotation and SIP distortion are all properties of the optics, not of which
    sub-region we kept.
    """
    out = header.copy()
    for key, shift in (("CRPIX1", x0), ("CRPIX2", y0)):
        if key in out:
            out[key] = float(out[key]) - shift
    return out


def _write(hdus, name):
    path = OUT / f"{name}.fits.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


def cut_cube(src, name, size, on_brightest=False):
    """Cut a native-layout stack (or a raw Bayer sub) and write it gzipped."""
    with fits.open(src) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
    if data.ndim == 3:
        centre = brightest_pixel(data) if on_brightest else None
        y0, y1, x0, x1 = centre_slice(data.shape[1:], size, centre)
        cut = data[:, y0:y1, x0:x1]
    else:  # raw Bayer: keep the cutout origin even so the mosaic phase is preserved
        y0, y1, x0, x1 = centre_slice(data.shape, size)
        y0 -= y0 % 2
        x0 -= x0 % 2
        y1, x1 = y0 + size, x0 + size
        cut = data[y0:y1, x0:x1]
    header = shift_wcs(header, x0, y0)
    header["SPCUTOUT"] = (f"[{x0}:{x1},{y0}:{y1}]", "cutout of the original frame")
    header["SPORIGIN"] = (Path(src).name, "source frame")
    path = _write([fits.PrimaryHDU(data=cut, header=header)], name)
    return path, (x0, y0, x1, y1), header


def cut_mef(src, name, size):
    """Cut a CrowdSky multi-extension frame, preserving its structure."""
    with fits.open(src) as hdul:
        primary = fits.PrimaryHDU(header=hdul[0].header.copy())
        planes, shape = [], None
        for hdu in hdul:
            ename = str(hdu.header.get("EXTNAME", "")).upper()
            if ename in ("RED", "GREEN", "BLUE") and hdu.data is not None:
                shape = hdu.data.shape
                planes.append((ename, hdu.data, hdu.header.copy()))
        star_tab = next((Table(h.data) for h in hdul
                         if str(h.header.get("EXTNAME", "")).upper() == "STAR-TAB"), None)
    y0, y1, x0, x1 = centre_slice(shape, size)
    primary.header = shift_wcs(primary.header, x0, y0)
    primary.header["SPCUTOUT"] = (f"[{x0}:{x1},{y0}:{y1}]", "cutout of the original frame")
    primary.header["SPORIGIN"] = (Path(src).name, "source frame")
    hdus = [primary]
    for ename, data, hdr in planes:
        hdus.append(fits.ImageHDU(data=data[y0:y1, x0:x1], header=hdr, name=ename))
    # Keep FOOTPRINT: it is the trap the loader must not mistake for a science plane, so
    # the bundled frame should still contain it.
    hdus.append(fits.ImageHDU(
        data=np.ones((y1 - y0, x1 - x0), dtype=np.uint8), name="FOOTPRINT"))
    if star_tab is not None:
        inside = ((np.asarray(star_tab["x"]) >= x0) & (np.asarray(star_tab["x"]) < x1)
                  & (np.asarray(star_tab["y"]) >= y0) & (np.asarray(star_tab["y"]) < y1))
        sub = star_tab[inside]
        sub["x"] = np.asarray(sub["x"]) - x0
        sub["y"] = np.asarray(sub["y"]) - y0
        hdus.append(fits.BinTableHDU(sub, name="STAR-TAB"))
    path = _write(hdus, name)
    return path, (x0, y0, x1, y1), primary.header


def write_sidecar(src_frame, name, x0, y0):
    """Write the cutout's ``.wcs`` sidecar from the original frame's cached solve.

    Bundling the solved WCS is what lets the docs and tests do real astrometry-dependent
    work -- cross-match, forced photometry, a light curve -- with no solver installed and
    no network.
    """
    sidecar = Path(src_frame).with_suffix(".wcs")
    if not sidecar.exists():
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        header = fits.getheader(sidecar)
    out = shift_wcs(header, x0, y0)
    path = OUT / f"{name}.wcs"
    fits.PrimaryHDU(header=out).writeto(path, overwrite=True)
    return path


def trim_gaia(footprints, size, margin_deg=0.05):
    """Gaia rows covering the bundled cutouts, with only the columns we need."""
    cat = Table.read(GAIA)
    keep = np.zeros(len(cat), dtype=bool)
    ra = np.asarray(cat["ra"], dtype=float)
    dec = np.asarray(cat["dec"], dtype=float)
    for header in footprints:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wcs = WCS(header)
        if not wcs.has_celestial:
            continue
        x, y = wcs.world_to_pixel_values(ra, dec)
        pad = margin_deg * 3600.0 / 2.39  # margin in pixels
        keep |= ((x > -pad) & (x < size + pad) & (y > -pad) & (y < size + pad))
    sub = cat[keep]
    cols = [c for c in GAIA_COLUMNS if c in sub.colnames]
    sub = sub[cols]
    sub.meta["comment"] = (
        "Gaia DR3 subset covering the bundled example cutouts, with synthetic "
        "Johnson-Kron-Cousins V. Trimmed from a 1.5 deg mosaic so the examples run offline."
    )
    path = OUT / "gaia_mwcam.ecsv"
    sub.write(path, format="ascii.ecsv", overwrite=True)
    return path, len(sub)


def pack(directory):
    """Tar+gzip the dataset and return ``(archive_path, sha256)``.

    Flat archive, sorted for reproducibility. The checksum is what ``examples.download``
    verifies, so a truncated or replaced asset fails loudly instead of surfacing later as
    inexplicable photometry.
    """
    import hashlib
    import tarfile

    directory = Path(directory)
    archive = directory / f"example-data-{DATA_VERSION}.tar.gz"
    files = sorted(p for p in directory.iterdir()
                   if p.is_file() and p.name not in (archive.name, "SHA256"))
    with tarfile.open(archive, "w:gz") as tf:
        for p in files:
            tf.add(p, arcname=p.name)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    (directory / "SHA256").write_text(
        f"{digest}  {archive.name}\n", encoding="utf-8")
    return archive, digest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--size", type=int, default=1000, help="cutout size in pixels")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    footprints, total = [], 0

    for name, pattern in STACKS.items():
        hits = sorted(glob.glob(str(S50 / pattern), recursive=True))
        if not hits:
            print(f"  MISSING {name}: {pattern}")
            continue
        src = hits[0]
        path, (x0, y0, _x1, _y1), header = cut_cube(
            src, name, args.size, on_brightest=(name == "stack_saturated"))
        side = write_sidecar(src, name, x0, y0)
        footprints.append(fits.getheader(side) if side else header)
        total += path.stat().st_size
        print(f"  {path.name:28s} {path.stat().st_size/1e6:5.2f} MB"
              f"  wcs={'yes' if side else 'NO'}")

    mef_hits = sorted(glob.glob(str(CROWDSKY / "*.fits")))
    if mef_hits:
        path, (x0, y0, _x1, _y1), header = cut_mef(mef_hits[0], "crowdsky_mef", args.size)
        footprints.append(header)
        total += path.stat().st_size
        # CrowdSky plate-solves server-side (PLTSOLVD), so its header WCS is
        # photometric-grade. Lift it into a sidecar as well, so `astrometry.load_wcs`
        # works uniformly across every bundled frame -- otherwise an example using the
        # MEF would try to plate-solve and demand an API key.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wcs = WCS(header, naxis=2)
        side = OUT / "crowdsky_mef.wcs"
        fits.PrimaryHDU(header=wcs.to_header(relax=True)).writeto(side, overwrite=True)
        total += side.stat().st_size
        print(f"  {path.name:28s} {path.stat().st_size/1e6:5.2f} MB  (MEF, WCS lifted)")

    raws = sorted(glob.glob(str(RAWS / "*.fit")))[:N_RAWS]
    for i, src in enumerate(raws, 1):
        path, _box, _h = cut_cube(src, f"raw_sub_{i}", args.size)
        total += path.stat().st_size
        print(f"  {path.name:28s} {path.stat().st_size/1e6:5.2f} MB  (raw Bayer sub)")

    gaia_path, n = trim_gaia(footprints, args.size)
    total += gaia_path.stat().st_size
    print(f"  {gaia_path.name:28s} {gaia_path.stat().st_size/1e6:5.2f} MB  ({n} sources)")

    print(f"\ntotal bundled: {total/1e6:.2f} MB in {OUT}")


if __name__ == "__main__":
    main()
