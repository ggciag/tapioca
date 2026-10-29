# How to install tapIOca

```bash
git clone https://github.com/ggciag/tapioca.git #or https://github.com/geojobueno/tapioca.git
cd tapioca
conda env create -f environment.yml -y
pip install -e .
```

`environment.yml` fixes the versions of the Python packages used to develop and test tapIOca, so a new
release of a library does not change the results.

The scenarios drawn in SVG files are read as **SVG 1.1** files saved by **Inkscape 1.4** (tested with
Inkscape 1.4.4). Inkscape is not installed by `environment.yml`: install it from
[inkscape.org](https://inkscape.org). If a file declares another version of SVG, or was saved by another
version of Inkscape, `tapioca-svg`, `tapioca-svg-interfaces` and `SVGScenario` print a warning (the versions are
`tapioca.inputs.svg.SVG_VERSION` and `INKSCAPE_VERSION`), and `infos_<folder>.txt` records the versions of the
file and of tapIOca used to build the scenario.

To check that everything works, run the tests (with pytest, included in `environment.yml`) from the root of the repository:

```bash
pytest
```