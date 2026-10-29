# Scenarios from SVG files

One example of each numerical experiment of `tapioca.inputs.svg.SVGScenario`. Open the SVGs in Inkscape
(*Object > Object Properties*) to see the label, title and description of each unit, and the description of
the `setup` rectangle (in the `domain` layer) with the domain, the resolution and the experiment. The units
are in the `lithologies` layer. Every unit gives its reference thickness (`h = ... km`), except the
asthenosphere, which fills the rest of the domain.

| Folder | Scenario | Experiment |
|---|---|---|
| `rifting` | continental lithosphere (155 km) with a weak seed | rifting, 1 cm/yr, magmatism on |
| `rifting_symmetric` | the same lithosphere, with the weak seed at the centre | rifting, 1 cm/yr, both boundaries moving outwards at 0.5 cm/yr (`moving_boundary = symmetric`) |
| `rifting_with_basin` | continental lithosphere thinned below a sedimentary basin, with a weak seed (lithologies drawn at half horizontal scale) | rifting, 1 cm/yr, magmatism on |
| `compression` | oceanic lithosphere - continental margin ramp - non-cratonic lithosphere - craton, with a dipping weak zone | compression, the oceanic plate pushed from the left at 1 cm/yr |
| `accordion` | non-cratonic lithosphere between two cratons, with a weak seed (M300) | accordion (extension, convergence, quiescence and extension) |
| `quiescence` | cratonic keel with a mobile belt in a 4000 km domain | quiescence, 1000 Myr |

To create the input files of Mandyoc:

```bash
cd rifting
tapioca-svg rifting.svg -o inputs
```

The settings of the SVG can be changed on the command line, e.g. `tapioca-svg accordion.svg --quiescence-time 0 --dx 1 --dz 1`.

In rifting and accordion, `moving_boundary = symmetric` (in the SVG, as in `rifting_symmetric`, or `--moving-boundary symmetric`)
puts the reference frame at the centre of the domain: both boundaries move outwards with half of the velocity,
e.g. `tapioca-svg rifting.svg -o inputs --moving-boundary symmetric`.
In compression, the opposite boundary of the moving one is fixed.

The `rifting` and `quiescence` SVGs have the line `preset_initial_temperature = ../preliminary` commented
out: uncomment it (or use `--preset-initial-temperature`) with the folder of a preliminary scenario to start
from its final thermal state instead of the steady-state conduction in the lithosphere.
