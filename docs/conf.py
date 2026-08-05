"""Sphinx configuration.

Markdown throughout (MyST), so the design-decision documents already written as ``.md``
are part of the site without conversion.
"""

project = "seestar-photometry"
author = "Kieran Leschinski"
copyright = "2026, Kieran Leschinski"

try:
    from importlib.metadata import version as _version

    release = _version("seestar-photometry")
except Exception:  # not installed (e.g. a docs-only checkout)
    release = "0.2.0"
version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",       # the codebase uses numpydoc-style docstrings
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
]

# "linkify" is deliberately absent: it auto-links bare URLs but needs
# linkify-it-py, which is not worth a docs dependency for that.
myst_enable_extensions = ["deflist", "colon_fence"]
myst_heading_anchors = 3

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "README.md"]

# Autodoc: document what the module docstrings say, in source order rather than
# alphabetically -- these modules are written to be read top to bottom.
# Render numpydoc "Attributes" sections as :ivar: rather than separate attribute
# directives. Without this, autodoc documents each dataclass field twice -- once from
# the class docstring and once from the annotation -- and Sphinx warns on every one.
napoleon_use_ivar = True

autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
# sep and astroalign are C extensions / optional extras; don't require them to build docs.
autodoc_mock_imports = ["sep", "astroalign", "skimage"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "astropy": ("https://docs.astropy.org/en/stable", None),
}

html_theme = "furo"
html_title = f"seestar-photometry {version}"
html_static_path = []
