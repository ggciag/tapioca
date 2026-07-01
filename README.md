# tapIOca

TapIOca is a python package with tools for pre- and post-processing of [Mandyoc](https://github.com/ggciag/mandyoc) models.

## About

The current version (v0.1.0) of tapIOca enables the manipulation of Mandyoc outputs from .nc files.
These netCDF files must be created from text file outputs containing the "x", "z" and "time" dimensions. Future versions of Mandyoc will natively support these output formats. In the meantime, you can use [Julia scripts](src/postrunning/README.md) to format your .txt data. You can send an e-mail to [jbueno@usp.br](mailto:jbueno@usp.br) in case of any question.

Future implementations of tapIOca are planned to include:
- A basic usage guide;
- Frameworks to facilitate scenarios creation;
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

## Contributing

This package is still in its early stages, so **contributions are welcome**. You can contribute by: sharing your own functions and methods, **reporting bugs or malfunctions**, suggesting new features and approaches for managing Mandyoc data.

## License

This is free software: you can redistribute it and/or modify it under the terms
of the **BSD 3-clause License**. A copy of this license is provided in
[LICENSE](https://github.com/ggciag/tapioca/blob/main/LICENSE).
