#!/usr/bin/env python
"""Cache the JKC and SDSS passband curves the S50 system is compared against.

Fetched once from the SVO Filter Profile Service and written to
``tools/standard_passbands.ecsv`` on the XP sampling grid, so
``build_s50_system_figures.py`` needs no network and the comparison is reproducible
offline. Transmissions are relative; only the shape matters, since every magnitude here
is an AB integral of the same spectra.

Johnson B/V and Cousins R/I rather than any other flavour of "JKC", because that is the
combination Gaia's ``synthetic_photometry_gspc`` reports and therefore what
``v_jkc_mag`` in :mod:`gaiadb` is on.

Usage::

    uv run python tools/fetch_standard_passbands.py

Maintenance script, not shipped in the wheel.
"""

import urllib.request
from pathlib import Path

import numpy as np

SVO = "http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id="

#: Column name -> SVO filter id. Order is blue to red within each system.
FILTERS = {
    "B_JKC": "Generic/Johnson.B",
    "V_JKC": "Generic/Johnson.V",
    "R_JKC": "Generic/Cousins.R",
    "I_JKC": "Generic/Cousins.I",
    "g_SDSS": "SLOAN/SDSS.g",
    "r_SDSS": "SLOAN/SDSS.r",
    "i_SDSS": "SLOAN/SDSS.i",
}

#: The XP product's grid, which everything is resampled onto: 336-1020 nm, 2 nm.
XP_GRID = np.arange(336.0, 1021.0, 2.0)

OUT = Path(__file__).resolve().parent / "standard_passbands.ecsv"


def fetch(filter_id):
    """``(wavelength_nm, transmission)`` for one SVO filter, peak-normalised."""
    with urllib.request.urlopen(SVO + filter_id, timeout=60) as response:
        text = response.read().decode()
    rows = np.array([[float(x) for x in line.split()]
                     for line in text.splitlines() if line.strip()
                     and not line.startswith("#")])
    wl_nm = rows[:, 0] / 10.0                      # SVO serves Angstrom
    return wl_nm, rows[:, 1] / rows[:, 1].max()


def main():
    from astropy.table import Table

    out = Table({"wavelength_nm": XP_GRID})
    for name, filter_id in FILTERS.items():
        wl, transmission = fetch(filter_id)
        # Zero outside the tabulated range: a filter that runs past the XP grid (I_JKC
        # does, slightly) is truncated rather than extrapolated, and that truncation is
        # stated in the doc rather than hidden.
        out[name] = np.interp(XP_GRID, wl, transmission, left=0.0, right=0.0)
        covered = (wl.min() <= XP_GRID[0]) or (wl.max() >= XP_GRID[-1])
        print(f"  {name:7s} {filter_id:22s} {wl.min():6.1f}-{wl.max():6.1f} nm"
              f"{'  (clipped by the XP grid)' if covered else ''}")

    out.meta["comments"] = [
        "JKC and SDSS passbands from the SVO Filter Profile Service, resampled onto",
        "the Gaia XP sampled grid (336-1020 nm, 2 nm) and peak-normalised.",
        "Johnson B/V + Cousins R/I: the combination Gaia GSPC reports.",
        f"Built by tools/{Path(__file__).name}.",
    ]
    out.write(OUT, overwrite=True)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
