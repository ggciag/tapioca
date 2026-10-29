"""Tests of the materials of the SVG scenarios (tapioca.inputs.svg.materials)."""

import numpy as np
import pytest

import tapioca as tp
from tapioca.inputs.svg import materials
from tapioca.inputs.svg.materials import RHEOLOGIES, LithologicalUnit, MandyocLayer, get_rheology


@pytest.mark.parametrize("key", list(tp.DEFAULT_MATERIALS))
def test_one_rheology_per_default_material(key):
    unit = RHEOLOGIES[key.lower()]
    assert unit.params == {**tp.DEFAULT_MATERIAL, **tp.DEFAULT_MATERIALS[key]}
    assert getattr(materials, key) is unit  # WET_OLIVINE, AIR, ... exported by the module


@pytest.mark.parametrize("name", ["wet olivine", "wet_olivine", "Wet Olivine", "WetOlivine", "WET_OLIVINE", "wet-olivine"])
def test_get_rheology_free_form_names(name):
    assert get_rheology(name) is materials.WET_OLIVINE


def test_get_rheology_unknown_name():
    with pytest.raises(ValueError, match="Unknown rheology"):
        get_rheology("granite")


def test_lithological_unit_defaults_and_flow_law():
    unit = LithologicalUnit("test", A=1.0e-20, n=3.0, Q=400e3, V=10e-6)
    assert unit.params["k"] == tp.DEFAULT_MATERIAL["k"]
    assert [unit.params[key] for key in "AnQV"] == [1.0e-20, 3.0, 400e3, 10e-6]


def test_lithological_unit_rejects_unknown_parameters():
    with pytest.raises(ValueError, match="Unknown parameter"):
        LithologicalUnit("test", viscosity=1.0e21)


def test_mandyoc_layer_material_overrides_the_rheology():
    layer = MandyocLayer("crust", "wet quartz", Nx=5, rho=2750.0, C=2.0, H=1.0e-10, cohesion_max=30e6)
    params = layer.material()
    assert layer.given == {"rho": 2750.0, "C": 2.0, "H": 1.0e-10, "cohesion_max": 30e6}
    assert set(params) == set(tp.DEFAULT_MATERIAL)
    assert (params["rho"], params["C"], params["H"], params["cohesion_max"]) == (2750.0, 2.0, 1.0e-10, 30e6)
    assert params["A"] == materials.WET_QUARTZ.params["A"]
    assert params["friction_angle_min"] == materials.WET_QUARTZ.params["friction_angle_min"]
    np.testing.assert_array_equal(layer.interface, np.zeros(5))


def test_mandyoc_layer_rejects_unknown_parameters():
    with pytest.raises(ValueError, match="Unknown parameter"):
        MandyocLayer("crust", "wet quartz", Nx=5, density=2750.0)  # descriptions accept the long names, layers do not


def test_mandyoc_layer_drawing_attributes():
    layer = MandyocLayer("crust", "wet quartz", Nx=5, title="Upper crust", color="#e49c7c", reference_thickness=25e3,
                         rho=2700.0)
    assert (layer.title, layer.color, layer.reference_thickness) == ("Upper crust", "#e49c7c", 25e3)
    assert MandyocLayer("crust", "wet quartz", Nx=5, rho=2700.0).title == "crust"  # default: the label
