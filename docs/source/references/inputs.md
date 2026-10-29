# Input builders

Input builders are in the subpackage `tapioca.inputs`. They create the input files of Mandyoc
(`interfaces.txt`, `param.txt`, `input_temperature_0.txt`, `input_velocity_0.txt`, ...). The generic
builders (`MandyocBuilder`, `VelocityFieldBuilder` and `TemperatureFieldBuilder`) are used from scripts,
and `tapioca.inputs.svg` builds a scenario drawn in an SVG file. See
[Scenarios from SVG files](../user_guide/svg_scenarios.md) for how to draw the scenario.

---

## Scenario builder

The `SVGScenario` class reads the SVG, builds the layers, the initial temperature and the velocity
boundary conditions of the numerical experiment (rifting, compression, accordion or quiescence), and writes
the input files and figures. The source code is in `inputs/svg/scenario.py`.

```{eval-rst}
.. currentmodule:: tapioca.inputs.svg.scenario

.. autosummary::
   :toctree: ../generated/

   SVGScenario
   SVGScenario.write_all
   SVGScenario.write_scripts
   SVGScenario.build_temperature
   SVGScenario.build_velocity
   SVGScenario.parameters
   SVGScenario.scale_bcv
   SVGScenario.infos
   SVGScenario.plot_setup
   SVGScenario.plot_temperature
   SVGScenario.yield_strength_envelope
```

Example:
```python
from tapioca.inputs.svg import SVGScenario

scenario = SVGScenario("rifting.svg", regime="accordion", quiescence_time=0)
scenario.write_all("my_scenario", particles_per_element=100)  # keyword arguments change param.txt
```

`SVGScenario` is built on the generic builders: `scenario.builder` is a `MandyocBuilder` with the interfaces
and the material fields of the SVG, the temperature is solved by a `TemperatureFieldBuilder` (`solve_heat_2D`, or
`solve_heat_1D` with `temperature_solver="columns"`) and the boundary velocities of the regime are set by
`VelocityFieldBuilder.set_kinematic_regime`. The settings are checked by `SVGChecks` (`inputs/checks.py`). The input
files are exported by `scenario.builder`, so they can be changed by script before writing:

```python
scenario = SVGScenario("rifting.svg")
scenario.build_temperature()
scenario.build_velocity()
scenario.velocity_builder.set_region("top", [700e3, 900e3], vz=1.0e-11, mode="add")  # m/s
scenario.builder.create_velocity_field(scenario.velocity_builder)
scenario.builder.calculate_YSE()            # 2D fields in scenario.builder.DTree.fields
scenario.write_all("my_scenario")
```

The same from the command line (installed with the package):
```bash
tapioca-svg rifting.svg -o my_scenario --regime accordion --quiescence-time 0
```

The run script of Mandyoc is written for the system that runs tapioca (`run_linux.sh` or `run_mac.sh`).
Another system, the number of MPI processes and the options of Mandyoc may be given:
```bash
tapioca-svg rifting.svg -o my_scenario --system linux --ncores 20 --mandyoc-options="-seed 0,5,8 -strain_seed 0.0,1.0,1.0"
```

---

## Geometry and materials from SVG files

Functions to read the setup shape and the units of the SVG and build the interfaces. The source code is in `inputs/svg/reader.py`.
`layers_to_builder` turns the layers of the SVG into a `MandyocBuilder`, and `write_interfaces` exports its `interfaces.txt`
(with the 13 parameters of `tapioca.DEFAULT_MATERIAL`).

```{eval-rst}
.. currentmodule:: tapioca.inputs.svg.reader

.. autosummary::
   :toctree: ../generated/

   read_setup
   read_units
   layers_from_svg
   layers_to_builder
   write_interfaces
   plot_layers
   layer_thickness
```

The errors of the SVG and of the settings are `SVGInputError` (`inputs/checks.py`, a `ValueError`).

---

## Preset of the initial temperature

Functions of `tapioca._aux_functions` used with `preset_initial_temperature`: the final thermal state of a
preliminary scenario (read with `MandyocScen`, or from its `temperature_*.txt`) adapted to the lithosphere of
each column of the new scenario.

```{eval-rst}
.. currentmodule:: tapioca._aux_functions

.. autosummary::
   :toctree: ../generated/

   read_preset_temperature
   adapt_preset_temperature
   read_interfaces
   convert_units
```

Only the interfaces from the command line:
```bash
tapioca-svg-interfaces scenario.svg -o interfaces.txt --plot layers.png
```

---

## Materials

Layers of a Mandyoc scenario and their rheologies. The source code is in `inputs/svg/materials.py`.
The rheologies are instances of `LithologicalUnit` built from `tapioca.DEFAULT_MATERIALS` (`WET_OLIVINE`,
`DRY_OLIVINE`, `WET_QUARTZ`, `WET_ANORTHITE`, `LINEAR_SALT`, `DRY_QUARTZ`, `BASALT`, `PLAGIOCLASE` and `AIR`);
a material added to that dictionary is also available in the SVG descriptions.

```python
from tapioca.inputs.svg import WET_OLIVINE, LithologicalUnit, get_rheology

WET_OLIVINE.params["V"]                       # parameters of interfaces.txt (rho, C, H, A, n, Q, V, k, ...)
get_rheology("wet olivine") is WET_OLIVINE    # free-form names, as in the SVG descriptions
my_olivine = LithologicalUnit("my_olivine", **{**WET_OLIVINE.params, "V": 15.0e-6})
```

```{eval-rst}
.. currentmodule:: tapioca.inputs.svg.materials

.. autosummary::
   :toctree: ../generated/

   MandyocLayer
   LithologicalUnit
   get_rheology
```
