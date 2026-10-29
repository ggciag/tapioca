"""Tests of the reader of SVG scenarios (tapioca.inputs.svg.reader)."""

import re

import numpy as np
import pytest

from tapioca.inputs.svg import INKSCAPE_VERSION, SVG_VERSION, SVGInputError, SVGScenario, svg_versions, layers_from_svg, layers_to_builder, read_setup, read_units, write_interfaces
from tapioca.inputs.svg.reader import eval_number, parse_description

from conftest import SVG_SCENARIOS, quiet, svg_path


def _layers(svg, Nx=161):
    setup = read_setup(str(svg))
    with quiet():
        return setup, layers_from_svg(str(svg), setup["Lx"], setup["Lz"], Nx, setup_id=setup["id"],
                                       min_thickness=setup.get("min_thickness"), verbose=False)


# Descriptions ##################################################################

def test_eval_number_accepts_only_arithmetic():
    assert eval_number("1.67e-6 / 2700.0") == pytest.approx(1.67e-6 / 2700.0)
    assert eval_number("2**3 - 1") == 7.0
    with pytest.raises(ValueError):
        eval_number("__import__('os')")


def test_parse_description():
    rheology, props, thickness = parse_description("wet quartz\nC = 1.0\nrho = 2700\nH = 1e-10\nh = 25 km\ncohesion_max = 3e7", "test")
    assert rheology.name == "wet_quartz"
    assert props == {"C": 1.0, "rho": 2700.0, "H": 1e-10, "cohesion_max": 3e7}
    assert thickness == 25e3


@pytest.mark.parametrize("text, message", [
    ("C = 1.0\nrho = 2700\nH = 0", "rheology"),
    ("wet quartz\nC = 1.0\nH = 0", "missing rho"),
    ("wet quartz\nC = 1.0\nrho = 2700\nH = 0\nviscosity = 1e21", "unknown property"),
    ("granite\nC = 1.0\nrho = 2700\nH = 0", "Unknown rheology"),
])
def test_parse_description_errors(text, message):
    with pytest.raises(SVGInputError, match=message):
        parse_description(text, "test")


# Examples ######################################################################

@pytest.mark.parametrize("name", SVG_SCENARIOS)
def test_example_setup(name):
    setup = read_setup(str(svg_path(name)))
    assert setup is not None, "the example must have a shape labelled 'setup'"
    assert setup["Lx"] > 0 and setup["Lz"] > 0 and setup["Nx"] > 1 and setup["Nz"] > 1
    regime, moving_boundary = SVG_SCENARIOS[name]
    assert setup["experiment"]["regime"] == regime
    assert setup["experiment"].get("moving_boundary", "right") == moving_boundary


@pytest.mark.parametrize("name", SVG_SCENARIOS)
def test_example_units_have_reference_thickness(name):
    '''Every unit of the examples gives h, except the bottom one (asthenosphere), which fills the rest of the domain.'''
    _, units = read_units(str(svg_path(name)))
    _, layers = _layers(svg_path(name))
    bottom = layers[0].layer_label
    without_h = [label for label, unit in units.items() if unit["thickness"] is None and label != bottom]
    assert without_h == []
    assert units[bottom]["thickness"] is None


@pytest.mark.parametrize("name", SVG_SCENARIOS)
def test_example_layers(name):
    setup, layers = _layers(svg_path(name))
    assert layers[0].layer_label == "asthenosphere"
    assert layers[-1].rheology_name == "air"
    bases = np.array([layer.interface for layer in layers])
    assert np.all(np.diff(bases, axis=0) <= 0), "the layers are stacked from the bottom to the top"
    assert np.all(bases[1:] >= 0) and np.all(bases[1:] <= setup["Lz"])


def test_h_sets_the_thickness_of_a_layer(tmp_path):
    '''Changing h of the upper crust of the rifting example changes its thickness, but not the drawing.'''
    text = svg_path("rifting").read_text()
    assert "h = 25 km" in text
    svg = tmp_path / "rifting.svg"
    svg.write_text(text.replace("h = 25 km", "h = 30 km", 1))
    _, layers = _layers(svg)
    labels = [layer.layer_label for layer in layers]
    i = labels.index("upper_crust")
    thickness = layers[i].interface - layers[i + 1].interface
    assert np.median(thickness) == pytest.approx(30e3)


def test_h_of_the_bottom_unit_is_an_error(tmp_path):
    text = svg_path("rifting").read_text()
    svg = tmp_path / "rifting.svg"
    svg.write_text(re.sub(r"(rho = 3378\.0\nH = 0\.0)", r"\1\nh = 100 km", text, count=1))
    with pytest.raises(SVGInputError, match="bottom layer"):
        _layers(svg)


# Builder and interfaces.txt ######################################################

def test_layers_to_builder(tmp_path):
    setup, layers = _layers(svg_path("rifting"))
    builder = layers_to_builder(layers, setup["Lx"], setup["Lz"], Nz=31, path=str(tmp_path))
    interfaces = builder.DTree.interfaces.ds
    assert len(interfaces.data_vars) == len(layers)
    ids = sorted(interfaces[name].attrs["id"] for name in interfaces.data_vars)
    assert ids == [-1] + list(range(len(layers) - 1))
    seed = [name for name in interfaces.data_vars if name.endswith("weak_seed")]
    assert interfaces[seed[0]].attrs["weakening_seed"] == 1.0  # initial strain in the weak seed
    others = [name for name in interfaces.data_vars if not name.endswith("weak_seed")]
    assert all(interfaces[name].attrs["weakening_seed"] == -1.0 for name in others)


def test_write_interfaces_does_not_touch_other_files(tmp_path):
    setup, layers = _layers(svg_path("rifting"))
    (tmp_path / "interfaces.txt").write_text("keep me")
    write_interfaces(layers, setup["Lx"], setup["Lz"], str(tmp_path / "other.txt"))
    assert (tmp_path / "interfaces.txt").read_text() == "keep me"
    lines = (tmp_path / "other.txt").read_text().splitlines()
    assert lines[0].split()[0] == "C" and len(lines) == 13 + 161


def test_only_the_setup_label_is_read(tmp_path):
    """The shape of the domain must be labelled setup (the older label domain is not read)."""
    svg = tmp_path / "old.svg"
    svg.write_text(svg_path("rifting").read_text().replace('inkscape:label="setup"', 'inkscape:label="domain"')
                   .replace('id="setup"', 'id="domain"'))
    assert read_setup(str(svg)) is None
    with pytest.raises(SVGInputError, match="No shape labelled 'setup'"):
        SVGScenario(str(svg))


@pytest.mark.parametrize("name", SVG_SCENARIOS)
def test_examples_are_svg_1_1(name):
    svg_version, inkscape_version = svg_versions(str(svg_path(name)))
    assert svg_version == SVG_VERSION
    assert inkscape_version is None or inkscape_version.startswith(INKSCAPE_VERSION + ".")


@pytest.mark.parametrize("old, new, warning", [
    ('version="1.1"', 'version="2.0"', "declares SVG 2.0"),
    ('<svg', '<svg inkscape:version="1.5.0 (abc, 2027-01-01)"', "saved by Inkscape 1.5.0"),
])
def test_other_versions_are_warned(tmp_path, capsys, old, new, warning):
    svg = tmp_path / "other.svg"
    svg.write_text(svg_path("rifting").read_text().replace(old, new, 1))
    scenario = SVGScenario(str(svg), dx=10, dz=10, verbose=False)
    assert warning in capsys.readouterr().out
    assert any("tested: SVG 1.1, Inkscape 1.4" in line for line in scenario.infos())


def test_names_of_interfaces_parameters_in_the_description():
    """The long names of tapioca.INTERFACES_PARAMETERS are the same parameters as their keys."""
    _, props, _ = parse_description("wet quartz\nscale_factor = 1.0\ndensity = 2700\nradiogenic_heat = 1e-10\n"
                                    "thermal_diffusivity = 1e-6\nweak_seed_strain = 0.5", "test")
    assert props == {"C": 1.0, "rho": 2700.0, "H": 1e-10, "k": 1e-6, "weakening_seed": 0.5}
