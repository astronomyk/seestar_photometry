"""Import-time contracts.

The lazy-matplotlib rule was documented for the predecessor package but never actually
tested, and had to be re-verified by hand. This is that test.
"""

import subprocess
import sys


def _run(code):
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def test_import_does_not_pull_in_matplotlib():
    """``import seestar_photometry`` must work in a core install with no plot extra.

    A fresh interpreter is used because the test session itself imports matplotlib for
    the report tests, which would make an in-process check pass regardless.
    """
    out = _run(
        "import sys, seestar_photometry; "
        "print('matplotlib' in sys.modules)"
    )
    assert out == "False"


def test_import_does_not_pull_in_astroquery():
    """astroquery is only needed to *build* a catalogue, not to use a cached one, and
    importing it is slow -- which matters when every worker process pays the cost."""
    out = _run(
        "import sys, seestar_photometry; "
        "print('astroquery' in sys.modules)"
    )
    assert out == "False"


def test_import_does_not_pull_in_the_catalog_extra():
    """``pyarrow`` and ``astropy_healpix`` back an opt-in multi-GB download.

    A core install has neither, and ``gaiadb`` is eagerly imported by the package, so
    every one of its imports has to sit inside a function. ``astroalign`` is here too:
    it was only the stacking extra until the local solver started bootstrapping on it.
    """
    out = _run(
        "import sys, seestar_photometry; "
        "print([m for m in ('pyarrow', 'astropy_healpix', 'astroalign') "
        "if m in sys.modules])"
    )
    assert out == "[]"


def test_plots_module_imports_without_matplotlib_at_module_level():
    """``plots`` is eagerly imported by the package, so it must stay import-light."""
    out = _run(
        "import sys; from seestar_photometry import plots, report; "
        "print('matplotlib' in sys.modules)"
    )
    assert out == "False"


def test_public_names_are_exported():
    import seestar_photometry as sp

    for name in sp.__all__:
        assert hasattr(sp, name), name


def test_bands_order_is_canonical():
    from seestar_photometry import BANDS

    assert BANDS == ("R", "G", "B")


def test_green_is_index_one():
    """Load-bearing convention: green is the science band at axis-0 index 1."""
    from seestar_photometry import BANDS

    assert BANDS.index("G") == 1
