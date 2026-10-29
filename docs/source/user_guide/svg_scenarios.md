# Scenarios from SVG files

The geometry of a Mandyoc scenario can be drawn in an SVG file (e.g. in [Inkscape](https://inkscape.org))
and turned into the input files of Mandyoc with `tapioca.inputs.svg`:

```bash
tapioca-svg scenario.svg -o my_scenario
```

or, in Python:

```python
from tapioca.inputs.svg import SVGScenario
SVGScenario("scenario.svg").write_all("my_scenario")
```

The SVG must be an SVG 1.1 file saved by Inkscape (tested with Inkscape 1.4.4): tapioca reads the object
properties of Inkscape (label, title and description) and the shapes with `svgelements`, so curves,
transforms and groups may be used. A file that declares another SVG version, or that was saved by another
version of Inkscape, is read with a warning: check the layers (e.g. with `tapioca-svg-interfaces --plot`), since
a new version could change how the drawing is stored. The versions of the file are recorded in
`infos_<folder>.txt`.

The output folder receives `interfaces.txt`, `param.txt`, `input_temperature_0.txt`,
`input_velocity_0.txt`, `vel_bc.txt`, `velz_bc.txt`, `scale_bcv.txt` (accordion), `infos_<folder>.txt`,
a copy of the SVG and the figures `numerical_setup.png` and `initial_temperature.png`.
Examples are in `examples/svg_scenarios`.

It also receives the scripts to run the scenario (`--no-scripts` skips them):

- `run_linux.sh` or `run_mac.sh` (`--system`, default: the system that runs tapioca): runs Mandyoc with
  `$PETSC_DIR/$PETSC_ARCH/bin/mpirun` and then moves the outputs to one folder per variable with
  `0_organize_outputs.sh`. On the Mac, `caffeinate` keeps the computer awake during the run. `PETSC_DIR`,
  `PETSC_ARCH` and `MANDYOC` (the executable) set in the environment override the values of the script;
  `--ncores` and `--mandyoc-options` set the number of MPI processes and the options of Mandyoc;
- `0_organize_outputs.sh`: copy of `src/postrunning/0_organize_outputs.sh`;
- `zipper.sh`: zips the inputs and the outputs in `<folder>.zip`, skipping the files already in the zip,
  so it can be run again to add new outputs.

## Units of rock

Each filled shape is a body of rock described by its object properties (Inkscape: *Object > Object Properties*):

- **Label**: name of the unit (`layer_label`). Shapes with the same label are the same unit.
- **Title**: name of the unit in the legend of the figures (optional).
- **Description**: the rheology in the first line, a material of `tapioca.DEFAULT_MATERIALS` (`wet olivine`,
  `dry olivine`, `wet quartz`, `wet anorthite`, `linear salt`, `dry quartz`, `basalt`, `plagioclase` or `air`)
  and the properties in the next lines:

```
wet quartz
C = 1.0
rho = 2800.0
H = 0.19e-6 / 2800.0
h = 10 km
```

`C` (viscosity scale factor), `rho` (density) and `H` (radiogenic heat production, W/kg) are required.
`k`, `weakening_seed`, `cohesion_min`, `cohesion_max`, `friction_angle_min` and `friction_angle_max` are also
accepted and, like `A`, `n`, `Q` and `V`, take the values of the rheology when they are not given:
`interfaces.txt` has the 13 parameters of `tapioca.DEFAULT_MATERIAL`. The keys are those of
`DEFAULT_MATERIAL` (`rho`, `C`, `H`, ...); their names in `tapioca.INTERFACES_PARAMETERS` (`density`,
`scale_factor`, `radiogenic_heat`, `thermal_diffusivity`, `weak_seed_strain`, ...) are also accepted.
`A`, `n`, `Q` and `V` are optional and override the values of the rheology in `DEFAULT_MATERIALS` (the
examples in `examples/svg_scenarios` give them, to keep the flow laws used before the materials were unified). A shape without description
inherits the one of the group (e.g. an Inkscape layer) that contains it; shapes without any description
(annotations, arrows) are ignored. The fill colour is used in the figures.

Shapes are painted in document order, as Inkscape renders them: a simple way to draw is a rectangle of
asthenosphere over the whole domain and the other units on top of it. Two labels are special:
`asthenosphere` (the bottom unit) and `weak_seed` (weak zones, which receive an initial strain of 1.0:
`weakening_seed = 1` in `interfaces.txt`, unless another value is given in the description).
The top unit must be the sticky air.

The shapes do not need to fill the setup shape exactly: points not covered by any unit (thin gaps between
shapes, shapes that stop short of the sides, the top or the bottom of the domain) take the unit of the
nearest covered point, first along the same depth and then, for depths without any unit, along the same
column. A warning gives the number of points filled and the largest distance, so that a forgotten shape
is not filled silently.

Shapes that overlap a little also leave thin layers: a crust drawn a bit below the base of the air leaves
a thin layer of the mantle painted under it, between the air and the crust. A layer of the stack thinner
than `min_thickness` (default 1 km) in every column is removed, and takes the unit right below it, when
its unit is thicker elsewhere. Units thin everywhere (e.g. a 1 km decollement) and the thin ends of thick
layers are kept. A warning lists the removed layers; set `min_thickness = 0 km` in the description of the
setup shape (or `--min-thickness 0`) to keep them.

## Thickness of the units (`h`)

The drawing does not need to be in vertical scale: `h` gives the reference thickness of a unit.

- **Layers** (units that reach both sides of the domain): their most common drawn thickness becomes `h`
  and the same factor is applied in every column, so thinned regions keep their drawn proportion.
- **Bodies** (seeds, basins, weak zones, keels, ...): without `h` they follow the scaling of the unit
  they are drawn over; with `h` they get that thickness and the unit around them absorbs the difference.
- The **bottom unit** (asthenosphere) fills the rest of the domain, so it can not have `h`.

## Domain and numerical experiment

A shape labelled `setup` (e.g. a rectangle without fill in its own Inkscape layer, named `domain` in the
examples, with the units in the layer `lithologies`) gives the model domain:
its bounding box is mapped to the domain and its description gives the dimensions, the resolution and,
optionally, the numerical experiment:

```
Lx = 1600 km
Lz = 400 km
dx = 2 km
dz = 2 km
regime = accordion
velocity = 1 cm/yr
moving_boundary = right
quiescence_time = 0 Myr
time_max = 500 Myr
magmatism = off
magmatism_extraction = off
magmatism_Ts0 = 1080
preset_initial_temperature = ../preliminary
temperature_solver = 2D
```

Only the shape labelled `setup` is read (older SVGs with a shape labelled `domain` must be renamed).

Every setting can also be given to `SVGScenario` (as a keyword argument, with or without units, e.g.
`SVGScenario("scenario.svg", quiescence_time="0 Myr")`) or on the command line (e.g. `--regime`, `--dx`),
which overrides the SVG. The units are those of `tapioca.UNITS` (km or m; Myr, kyr or yr; cm/yr or mm/yr; oC).
Unknown settings and invalid values are errors (`SVGInputError`). Defaults:

| Setting | Default | |
|---|---|---|
| `regime` | `rifting` | `rifting` (extension), `compression` (convergence), `accordion` or `quiescence` |
| `velocity` | 1 cm/yr | velocity of the moving boundary (with `symmetric`, the sum of both boundaries) |
| `moving_boundary` | `right` | `left`, `right` or `symmetric` (rifting and accordion) |
| `quiescence_time` | 0 Myr | accordion and quiescence regimes |
| `time_max` | 60 Myr | rifting and compression; accordion: end of the second extension; quiescence: `quiescence_time` (or 1000 Myr) |
| `magmatism`, `magmatism_extraction` | `off` | with both `on`, a basalt layer receives the extracted melt (see below) |
| `magmatism_Ts0` | 1080 °C | |
| `magmatic_layer` | basalt at the top | label (or ID from the bottom) of the layer that receives the extracted melt |
| `water_layer` | none | label (or ID) of the unit drawn as water, right below the sticky air |
| `preset_initial_temperature` | none | folder of a preliminary scenario (see below) |
| `temperature_solver` | `2D` | `2D` or `columns` (see below) |
| `min_thickness` | 1 km | thin layers left by overlapping shapes (see above); `0 km` keeps them |

The accordion regime (as in the M300 experiments) scales the boundary velocity in time (`scale_bcv.txt`):
extension until 25 Myr, convergence until 75 Myr, quiescence during `quiescence_time` and extension again.

In the rifting and accordion regimes, the velocity is constant from the surface down to 20 km below the
thickest lithosphere of the domain; in the compression regime, down to 20 km below the lithosphere at the
moving boundary (e.g. an oceanic plate pushed towards a continent), and the opposite boundary is fixed.

With `moving_boundary = symmetric` (rifting and accordion), the reference frame is at the centre of the
domain: both boundaries move outwards with `velocity`/2, each with its own return flow below the constant
velocity, so that the integral of the inflow is zero on each boundary. The profiles are set by
`VelocityFieldBuilder.set_kinematic_regime`.

## Water and magmatic layers

tapioca does not guess which layer is the water or which one receives the melt: they are given by label (or by
the ID of the layer, from the bottom) in `water_layer` and `magmatic_layer`.

- `water_layer`: a unit drawn as water (e.g. with `air` rheology and `rho = 1000.0`), right below the sticky air.
  The surface (top of the solid layers), used by the temperature and the velocity, is then its base.
- `magmatic_layer`: by default, with `magmatism` and `magmatism_extraction` on, a basalt layer of zero thickness
  is added at the top of the solid layers: below the water when there is water, otherwise below the sticky air.
  A layer given in `magmatic_layer` (e.g. a basalt drawn in the SVG) is used as it is, and no basalt layer is
  added.

## Initial temperature

Without `preset_initial_temperature`, the initial temperature is the steady-state conduction in the lithosphere
(`TemperatureFieldBuilder`): from 0 °C at the surface to the mantle potential temperature (1350 °C) at the top of
the asthenosphere (LAB), with the conductivity and the radiogenic heat of each unit, and the mantle adiabat below
the LAB. With `temperature_solver = 2D` (default) the whole lithosphere is solved at once (`solve_heat_2D`), so heat
also flows sideways where the lithosphere changes along x (e.g. at the edge of a craton); with
`temperature_solver = columns` each column is solved alone (`solve_heat_1D`).

With `preset_initial_temperature` (folder with `temperature.nc`, read with `MandyocScen`, or with `param.txt` and
`temperature_*.txt`, or a NetCDF file with the temperature), the initial temperature is the horizontal mean of the final temperature of that
preliminary scenario, adapted to the lithosphere of each column: the columns with the lithosphere of the
preliminary scenario (read from its `interfaces.txt`) keep the profile; the others (keels, thinned or oceanic
lithosphere) keep it in the crust, get a linear geotherm in the lithospheric mantle down to the base of the
lithosphere and the asthenospheric profile below it.
