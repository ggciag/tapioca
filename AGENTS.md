# tapioca

Python package (`tapioca-mandyoc`) to build inputs for and post-process outputs of the Mandyoc geodynamic code, based on xarray and DataTrees.

## Layout

- `src/tapioca/scenClasses.py`: `MandyocScen`, which reads the Mandyoc outputs (`.nc` files).
- `src/tapioca/_variables.py`: constants (`SEC_PER_YEAR`, `UNITS`), `DEFAULT_MATERIAL`/`DEFAULT_MATERIALS`, `PARAMS_TYPES`.
- `src/tapioca/_aux_functions.py`: helpers (params, data and interfaces readers, units, heat solvers, yield strength, preset temperature).
- `src/tapioca/post_processing.py`, `src/tapioca/mandyocPlotter.py`: post-processing and plotting.
- `src/tapioca/inputs/`: the generic part, so there is no `inputs.generic`.
  - `generic_input_classes.py`: `GeometryDomain` and `VerboseExtension` (base classes), `MandyocBuilder`, `TemperatureFieldBuilder`, `VelocityFieldBuilder`.
  - `checks.py`: `InputChecks` and `SVGChecks` (checks of the input parameters), `InputError`, `SVGInputError`.
- `src/tapioca/inputs/svg/`: SVG scenarios (Inkscape drawings): `reader`, `materials`, `scenario` (`SVGScenario`), `cli`, `_variables` (regimes, default `param.txt` values, run scripts).
  - SVG scenarios are built on the generic builders: `SVGScenario.builder` is a `MandyocBuilder`, and every input file is exported through it.
  - Code specific to SVG belongs here, not in the generic builders. Generic helpers go to `_aux_functions.py`.
- `src/postrunning/`: scripts to organize and convert the outputs of Mandyoc, some copied next to the scenario inputs.
- `examples/svg_scenarios/`: example SVGs. Every layer has a reference thickness `h`, except the bottom one (the asthenosphere), which fills the rest of the domain.
- `docs/source/abstractions/svg_scenarios.md`: architecture of the SVG scenarios (Mermaid diagrams).
- `tests/`: pytest suite.

## Conventions

- Mandyoc grid: z is positive upwards (from `-Lz` at the base to `0` at the top) and node (i, j) = (0, 0) is the bottom-left corner. The SVG reader works with depths positive downwards and converts them in `layers_to_builder`.
- `MandyocBuilder` interfaces: ID `-1` is the top layer (sticky air). Interface `i` is the top of layer `i` and carries its material.
- `interfaces.txt` holds the 13 material parameters (`MATERIAL_PARAMETERS`).
- The weak seed layer has `weakening_seed = 1`; every other layer has `-1`.
- The water and magmatic layers are chosen by the user (`water_layer`, `magmatic_layer`, by label or ID), never guessed from their position.
- Initial temperature: steady-state conduction in the lithosphere (`solve_heat_2D` by default, `solve_heat_1D` with `temperature_solver = columns`), with the LAB at 1350 oC.
- `SEC_PER_YEAR` is 365.25 days for every conversion.
- SVG files: SVG 1.1 saved by Inkscape 1.4 (`reader.SVG_VERSION`, `reader.INKSCAPE_VERSION`); other versions are read with a warning.
- The default materials come from `DEFAULT_MATERIALS`. The example SVGs write their own flow-law values (A, n, Q, V) in the layer descriptions. `DRY_QUARTZ` and `PLAGIOCLASE` have zero flow laws on purpose until their values are chosen.
- Code, docstrings, commits and docs are written in English.

## Working with other branches

Several branches (for example `particles`, `inputs`, `julia_in_chunks`) work on the shared core files at the same time:
`_aux_functions.py`, `_variables.py`, `mandyocPlotter.py`, `post_processing.py`, `scenClasses.py`, `inputs/generic_input_classes.py`.

- Keep edits to these files minimal. Prefer adding functions or methods over rewriting existing ones.
- Avoid renames, reformatting and cosmetic edits in these files.
- Bugs found outside the SVG code are flagged with `## BUG (Claude)` / `## NOTE (Claude)` comments for the owner of the generic classes to review; a fix that the SVG code needs is marked `## FIX (Claude)`.
- Before changing them, run `git fetch origin` and check what other branches changed.
- Before a commit or PR, simulate the merges, for example `git merge-tree --write-tree HEAD origin/<branch>`, and report conflicts.
- A fix to code that another branch is rewriting may belong in that branch instead.

## Commands

```bash
conda env create -f environment.yml   # environment "tapioca", with fixed versions
pip install -e ".[test]"
pytest                                 # about 10 s
```

Mandyoc runs are done on the group's cluster, Hypatia.
