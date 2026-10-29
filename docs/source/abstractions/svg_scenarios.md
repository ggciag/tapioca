# SVG Scenarios Architecture

`tapioca.inputs.svg` turns a scenario drawn in an SVG file (SVG 1.1, saved by Inkscape) into the input files of
Mandyoc. It reads the drawing and the settings, and leaves the physics to the generic builders of
`tapioca.inputs`: the SVG code builds the stack of layers, and `MandyocBuilder`, `TemperatureFieldBuilder` and
`VelocityFieldBuilder` build the fields and export the files.

## 1. Class Architecture

`SVGScenario` extends the mixins of the generic inputs (`SVGChecks` for the checks of the parameters,
`VerboseExtension` for the prints) and composes the generic builders. The reader is a module of functions that
returns the stack of layers (`MandyocLayer`), whose materials come from `LithologicalUnit`.

```mermaid
classDiagram
    %% Mixins of tapioca.inputs
    class VerboseExtension {
        <<Mixin>>
        +bool verbose
        +_print_verbose(message)
    }
    class InputChecks {
        <<Mixin>>
        +error = InputError
        +_check_choice()
        +_check_on_off()
        +_check_number()
        +_check_path()
        +_check_keys()
        +_check_domain()
    }
    class SVGChecks {
        +error = SVGInputError
        +_check_setup()
        +_check_layers()
        +_check_experiment()
    }
    InputChecks <|-- SVGChecks

    %% SVG scenario
    class SVGScenario {
        +list layers
        +MandyocBuilder builder
        +TemperatureFieldBuilder temperature_builder
        +VelocityFieldBuilder velocity_builder
        +build_temperature()
        +build_velocity()
        +parameters()
        +write_all()
        +plot_setup()
    }
    SVGChecks <|-- SVGScenario
    VerboseExtension <|-- SVGScenario

    %% Reader (module of functions)
    class reader {
        <<module>>
        +read_setup(svg)
        +read_units(svg)
        +layers_from_svg(svg, Lx, Lz, Nx)
        +layers_to_builder(layers, Lx, Lz, Nz)
        +write_interfaces(layers, Lx, Lz, fname)
    }
    class MandyocLayer {
        +str layer_label
        +LithologicalUnit rheology
        +ndarray interface
        +material() dict
    }
    class LithologicalUnit {
        +str name
        +dict params
    }
    reader ..> MandyocLayer : creates
    MandyocLayer --> LithologicalUnit
    SVGScenario ..> reader : uses

    %% Generic builders
    class MandyocBuilder {
        +DataTree DTree
        +create_interface()
        +create_materials_fields()
        +set_params()
        +export_interfaces()
        +export_field()
    }
    class TemperatureFieldBuilder {
        +DataArray temperature
        +apply_linear_temperature()
        +solve_heat_1D()
        +solve_heat_2D()
    }
    class VelocityFieldBuilder {
        +DataTree velocs
        +set_kinematic_regime()
        +conservate_veloc()
    }
    SVGScenario *-- MandyocBuilder
    SVGScenario *-- TemperatureFieldBuilder
    SVGScenario *-- VelocityFieldBuilder
    TemperatureFieldBuilder --> MandyocBuilder : fields
    VelocityFieldBuilder --> MandyocBuilder : fields
```

### Class Responsibilities

*   **`SVGChecks` (extends `InputChecks`, `inputs/checks.py`):** checks the settings and the drawing (setup shape,
    stack of layers, combinations of settings) and raises `SVGInputError` with a message that names the parameter.
*   **`SVGScenario` (`inputs/svg/scenario.py`):** reads the settings (argument > SVG > default), chooses the water and
    magmatic layers, and drives the generic builders: the temperature (`solve_heat_2D` by default, or `solve_heat_1D`),
    the velocity of the kinematic regime (`set_kinematic_regime`) and the files (`write_all`).
*   **`reader` (`inputs/svg/reader.py`):** turns the shapes into a single stack of layers, bottom to top, with the
    depth of the base of each layer in every column, and converts it into the interfaces of a `MandyocBuilder`.
*   **`MandyocLayer` / `LithologicalUnit` (`inputs/svg/materials.py`):** a layer of the stack and its material (the 13
    parameters of `DEFAULT_MATERIAL`, from `DEFAULT_MATERIALS` and the description of the shape).
*   **`MandyocBuilder`, `TemperatureFieldBuilder`, `VelocityFieldBuilder` (`inputs/generic_input_classes.py`):** the
    generic builders, which own the `DataTree` with the interfaces and the fields and export the input files.

---

## 2. From the Drawing to the Input Files

```mermaid
flowchart TD
    SVG["scenario.svg<br/>(SVG 1.1, Inkscape)"]

    subgraph reader["reader.layers_from_svg"]
        SETUP["read_setup<br/>shape labelled setup:<br/>Lx, Lz, dx, dz, settings"]
        UNITS["read_units<br/>label, title, description<br/>of each shape"]
        RASTER["rasterize<br/>unit of each sample point"]
        GAPS["fill_gaps, remove_slivers"]
        STACK["build_stack<br/>one stack for all columns"]
        SCALE["rescale_thickness<br/>h of each unit"]
        LAYERS["MandyocLayer list<br/>(bottom to top, depth of the base)"]
    end

    subgraph scenario["SVGScenario"]
        CHECKS["SVGChecks"]
        EXTRA["water_layer, magmatic_layer<br/>(basalt of zero thickness)"]
        TOB["layers_to_builder"]
    end

    MB["MandyocBuilder<br/>/interfaces, /fields"]
    TB["TemperatureFieldBuilder"]
    VB["VelocityFieldBuilder"]
    FILES["interfaces.txt, param.txt,<br/>input_temperature_0.txt, input_velocity_0.txt,<br/>vel_bc.txt, velz_bc.txt, scale_bcv.txt"]

    SVG --> SETUP --> CHECKS
    SVG --> UNITS --> RASTER --> GAPS --> STACK --> SCALE --> LAYERS
    LAYERS --> CHECKS --> EXTRA --> TOB --> MB
    MB --> TB --> MB
    MB --> VB --> MB
    MB --> FILES
```

### Steps

*   **Shapes to sample points:** `read_units` reads the object properties of each shape (Inkscape: *Object > Object
    Properties*): the first line of the description is the rheology, the next lines are `key = value` (`C`, `rho`
    and `H` are required; `h` is the reference thickness). `rasterize` maps the bounding box of the setup shape to the
    domain and samples the drawing (`Nx` columns, rows of 100 m), painting the shapes in document order as Inkscape
    does; it also keeps the unit painted over (`under`), the host of seeds and basins. `fill_gaps` fills the points
    not covered by any shape and `remove_slivers` removes thin layers left by shapes that overlap a little.
*   **Columns to one stack:** each interface of Mandyoc is a single line `z(x)`, so one stack must fit every column.
    `build_stack` joins the sequences of units of the columns with the shortest common supersequence (`scs`), e.g.
    `[asthenosphere, mantle, crust, air]` and `[asthenosphere, mantle, crust, seed, crust, air]` give the second one,
    with zero thickness for the seed and the upper crust slot where they are absent. Away from a body, its
    interfaces stay flat, so Mandyoc draws a box and not a funnel.
*   **Thickness `h`:** `rescale_thickness` makes the drawing independent of the vertical scale: a unit with `h` is
    stretched by `h / (its most common drawn thickness)`, bodies without `h` follow their host, and bodies with `h`
    get that thickness while their host absorbs the difference. The bottom unit fills the rest of the domain.
*   **Stack to builders:** `layers_to_builder` creates the interfaces of `MandyocBuilder` (see below) and
    `create_materials_fields` fills the material fields; the temperature and velocity builders add their fields.

---

## 3. Conventions

The reader works with depths positive downwards, as the SVG; the generic builders use the convention of Mandyoc:
z positive upwards, from `-Lz` at the base to `0` at the top, and node `(0, 0)` at the bottom-left corner.

```mermaid
flowchart LR
    L["layers[i].interface<br/>depth of the base of layer i (m, down)"] -->|"z = -depth"| I["interface ID i = top of layer i<br/>z = -layers[i+1].interface<br/>ID -1: sticky air"]
    I --> F["/fields: rho, C, H, A, n, Q, V, k, ...<br/>(z, x), row 0 at the base"]
```

*   **Interfaces:** ID `i` is the top of layer `i` and carries its material; ID `-1` is the sticky air.
*   **Surface and LAB:** the surface is the top of the solid layers (the base of the sticky air, or of the water
    when `water_layer` is given); the LAB is the top of the asthenosphere. Both are passed by ID to
    `TemperatureFieldBuilder` and to `set_kinematic_regime`.
*   **Water and magmatic layers:** chosen by the user (label or ID), never guessed from their position.
