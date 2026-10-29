"""Tests of the variables (tapioca._variables) and of the auxiliary functions (tapioca._aux_functions)."""

import numpy as np
import pytest

import tapioca as tp
from tapioca._aux_functions import convert_units, evaluate_params, read_data, read_interfaces, read_params, spsolve_heat_col

from conftest import CONTINENTAL_RIFT


# Variables ##################################################################

def test_material_parameters_have_units_and_names():
    assert tp.MATERIAL_PARAMETERS == list(tp.DEFAULT_MATERIAL)
    for key in tp.MATERIAL_PARAMETERS:
        assert key in tp.PARAMETERS_UNITS
        assert key in tp.INTERFACES_PARAMETERS


@pytest.mark.parametrize("name", list(tp.DEFAULT_MATERIALS))
def test_default_materials_have_every_parameter(name):
    assert set(tp.DEFAULT_MATERIALS[name]) == set(tp.DEFAULT_MATERIAL)


def test_seconds_per_year():
    assert tp.SEC_PER_YEAR == pytest.approx(365.25 * 24 * 3600)


# param.txt ###################################################################

def test_read_params_of_example():
    params = read_params(str(CONTINENTAL_RIFT / "param.txt"))
    assert params["nx"] == 121 and params["nz"] == 31
    assert params["lx"] == pytest.approx(1.2e6)
    for key, value in params.items():
        assert isinstance(value, tp.PARAMS_TYPES[key])


def test_evaluate_params_from_text_and_keywords():
    text = """
    # comment line
    nx = 11        # inline comment
    time_max = 1.0e6
    solver = direct
    """
    assert evaluate_params(params_str=text) == {"nx": 11, "time_max": 1.0e6, "solver": "direct"}
    assert evaluate_params(nx="11", magmatism_Ts0="1080") == {"nx": 11, "magmatism_Ts0": 1080.0}


def test_evaluate_params_rejects_unknown_parameters():
    with pytest.raises(KeyError):
        evaluate_params(not_a_parameter=1)


# Files of Mandyoc ##############################################################

def test_read_data_of_a_step(tmp_path):
    """Mandyoc writes the nodes with x varying fastest, after two header lines."""
    Nx, Nz = 4, 3
    field = np.arange(Nx * Nz, dtype=float).reshape(Nz, Nx)  # (Nz, Nx), row 0 is the base
    (tmp_path / "temperature_10.txt").write_text("# header\n# header\n" + "\n".join(map(str, field.ravel())) + "\n")
    np.testing.assert_array_equal(read_data(str(tmp_path / "temperature_10.txt"), Nx, Nz), field)


def test_read_interfaces_of_example():
    params, interfaces = read_interfaces(CONTINENTAL_RIFT / "interfaces.txt")
    assert list(params) == ["C", "rho", "H", "A", "n", "Q", "V"]
    assert all(len(v) == 7 for v in params.values())  # one value per layer
    assert interfaces.shape == (121, 6) and np.all(interfaces <= 0)  # z of the interfaces, negative downwards


@pytest.mark.parametrize("text, quantity, value", [("120 km", "length", 120e3), ("120e3 m", "length", 120e3),
                                                    ("300 Myr", "time", 300e6), ("5 mm/yr", "velocity", 0.5),
                                                    ("1080 oC", "temperature", 1080.0)])
def test_convert_units(text, quantity, value):
    assert convert_units(text, quantity) == pytest.approx(value)


def test_convert_units_needs_a_unit():
    with pytest.raises(ValueError, match="unit"):
        convert_units("120", "length")
    assert convert_units("120", "length", default="km") == 120e3


# Steady-state conduction ######################################################

def _column(H=0.0, k=3.0, N=121, thickness=120e3):
    """Index 0 is the base of the lithosphere (LAB) and index -1 the surface, as in Mandyoc (z upwards)."""
    dz = thickness / (N - 1)
    return dz, spsolve_heat_col(np.zeros(N), np.full(N, k), np.full(N, H), 0.0, 1400.0, dz)


def test_steady_column_without_heat_is_linear():
    dz, T = _column()
    assert T[0] == pytest.approx(1400.0) and T[-1] == pytest.approx(0.0)
    np.testing.assert_allclose(T, np.linspace(1400.0, 0.0, len(T)), atol=1e-8)


def test_steady_column_with_radiogenic_heat():
    """Analytic solution with constant k and H: T(d) = Tlab d / L + H d (L - d) / (2 k), d = depth."""
    H, k, L = 1.0e-6, 3.0, 120e3
    dz, T = _column(H=H, k=k, thickness=L)
    d = np.linspace(L, 0.0, len(T))
    np.testing.assert_allclose(T, 1400.0 * d / L + H * d * (L - d) / (2 * k), rtol=1e-8)
