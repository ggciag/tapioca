# tapIOca

TapIOca is a python package with tools for pre- and post-processing of [Mandyoc](https://github.com/ggciag/mandyoc) models.

## About

The current version (v0.1.0) of tapIOca enables the manipulation of Mandyoc outputs from .nc files.
These netCDF files must be created from text file outputs containing the "x", "z" and "time" dimensions. Future versions of Mandyoc will natively support these output formats. In the meantime, you can use [Julia scripts](src/postrunning/README.md) to format your .txt data. You can send an e-mail to [jbueno@usp.br](mailto:jbueno@usp.br) in case of any question.

TapIOca also creates the input files of Mandyoc from scenarios drawn in SVG files (e.g. in Inkscape), with the
geometry, the materials and the numerical experiment (rifting, compression, accordion or quiescence):

```
tapioca-svg scenario.svg -o my_scenario
```

See [examples/svg_scenarios](examples/svg_scenarios/README.md) and the user guide (*Scenarios from SVG files*).
The SVG files must be SVG 1.1 saved by Inkscape 1.4 (tested with Inkscape 1.4.4); other versions are read with a
warning, since a new version could change how the drawing is stored.

Future implementations of tapIOca are planned to include:
- A basic usage guide;
- Frameworks to facilitate scenarios creation (generic builders in `tapioca.inputs` and scenarios from SVG files in `tapioca.inputs.svg`);
- Specialised plotting functions (xarray and matplotlib functions are working);
- New post-processing tools;

## How to install

The current recommended way to install tapIOca is cloning this repository, creating a conda environment and installing the source code:

```
git clone https://github.com/ggciag/tapioca.git
cd tapioca
conda env create -f environment.yml -y
pip install -e .
```

## Tests

The tests use [pytest](https://docs.pytest.org) and cover the variables and auxiliary functions, the materials,
the generic builders (`tapioca.inputs`: temperature solvers and kinematic regimes), the checks of the input
parameters, the SVG scenarios of `examples/svg_scenarios` (on a coarse grid, with water and magmatic layers and
the preset of the initial temperature), the command line tools and the post-processing of
`examples/continental_rift`. From the root of the repository:

```
pip install -e ".[test]"   # or: conda install pytest
pytest
```

Run them before a pull request, and add tests for new functions in the `tests` folder.

## Contributing

This package is still in its early stages, so **contributions are welcome**. You can contribute by: sharing your own functions and methods, **reporting bugs or malfunctions**, suggesting new features and approaches for managing Mandyoc data.

## License

This is free software: you can redistribute it and/or modify it under the terms
of the **BSD 3-clause License**. A copy of this license is provided in
[LICENSE](https://github.com/ggciag/tapioca/blob/main/LICENSE).
