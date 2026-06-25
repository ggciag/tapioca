# Tapioca

**TapIOca is a python package with tools for pre- and post-processing for [Mandyoc](https://github.com/ggciag/mandyoc) models.**

## About

The current version (v0.1.0) of tapIOca enables the manipulation of Mandyoc outputs from .nc files.
These netCDF files must be created from text file outputs containing the "x", "z" and "time" dimensions. Future versions of Mandyoc will natively support these output formats. In the meantime, you can use [Julia scripts](src/postrunning/README.md) to format your .txt data. You can send an e-mail to [jbueno@usp.br](mailto:jbueno@usp.br) in case of any question.

Future implementations of tapIOca are planned to include:

* Expansion and improvement of the documentation;
    * A basic usage guide
    * Make Abstractions section clear
    * Improving/changing the documentation theme (API Refence is confuse) 
* Frameworks to facilitate scenarios creation;
    * Read a scenario as input
    * Create classes to facilitate scenario creation
    * Create structures to abstract materials and interfaces
    * *Read a drawing (such as an SVG file) to create interfaces*
* Specialised plotting functions and classes;
    * Classes to plot scenarios (artists, plotters, etc)
    * Functions to plot initial conditions and input parameters
* New post-processing tools;
    * Heat flux calculation
    * Improve $\tau$ function



## How to contribute

This package is still in its early stages, so **contributions are welcome**. You can contribute by: sharing your own functions and methods, **reporting bugs or malfunctions**, suggesting new features and approaches for managing Mandyoc data.

Feel free to open new issues at [https://github.com/ggciag/tapioca/issues](https://github.com/ggciag/mandyoc/issues).

## License

This is free software: you can redistribute it and/or modify it under the terms
of the **BSD 3-clause License**. A copy of this license is provided in
[LICENSE](https://github.com/ggciag/tapioca/blob/main/LICENSE).