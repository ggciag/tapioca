# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# Path setup
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('../../src'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'TapIOca'
author = 'João Bueno'
release = 'v0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",           # Allows Markdown parsing
    "sphinx.ext.autodoc",    # Pulls docstrings from Python modules
    "sphinx.ext.napoleon",   # Understands the NumPy docstring format
    "sphinx.ext.autosummary", # Creates the neat Xarray-style API tables
    "sphinx.ext.intersphinx", # Link objects and classes to external docs
    'sphinx.ext.mathjax'
]

# Parameters
napoleon_numpy_docstring = True  
napoleon_use_param = False       # Creates a cleaner, table-like look for parameters
napoleon_use_ivar = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/api/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}

myst_enable_extensions = ["dollarmath", "amsmath"]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = []

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
#html_static_path = ['_static']


"""
myst_enable_extensions = [
    "alert",
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "gfm_autolink",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]
"""