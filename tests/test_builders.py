"""Tests of the generic builders (tapioca.inputs: MandyocBuilder, VelocityFieldBuilder and TemperatureFieldBuilder)."""

import os

import numpy as np
import pytest

import tapioca as tp
from tapioca.inputs import MandyocBuilder, TemperatureFieldBuilder, VelocityFieldBuilder
from tapioca.inputs.svg import AIR, DRY_OLIVINE, WET_OLIVINE, WET_QUARTZ

from conftest import quiet

Lx, Lz, Nx, Nz = 1000e3, 300e3, 101, 31  # 10 km x 10 km


@pytest.fixture
def builder(tmp_path):
    '''
    Air (40 km), crust (35 km), lithospheric mantle (thinner on the right) and asthenosphere.
    Each interface is the top of the layer whose material it carries.
    '''
    b = MandyocBuilder(str(tmp_path / "scenario"), Nx=Nx, Nz=Nz, Lx=Lx, Lz=Lz, verbose=False)
    b.create_interface(-1, 0.0, "air", **AIR.params)
    b.create_interface(0, np.where(b.x < Lx / 2, -160e3, -120e3), "asthenosphere", **WET_OLIVINE.params)
    b.create_interface(1, -75e3, "mantle", **DRY_OLIVINE.params)
    b.create_interface(2, -40e3, "crust", **WET_QUARTZ.params)
    b.create_materials_fields()
    return b


# MandyocBuilder ################################################################

def test_builder_grid(builder):
    assert builder.x[0] == 0 and builder.x[-1] == Lx and len(builder.x) == Nx
    assert builder.z[0] == -Lz and builder.z[-1] == 0 and len(builder.z) == Nz


def test_export_interfaces_writes_in_the_builder_folder(builder):
    builder.export_interfaces(export="none")
    assert os.path.isfile(os.path.join(builder.path, "interfaces.txt"))


def test_interface_needs_a_density(builder):
    with pytest.raises(ValueError, match="density"):
        builder.create_interface(3, -50e3, "no_density", C=1.0)


def test_material_fields_follow_the_interfaces(builder):
    rho = builder.DTree.fields.rho
    left, right = rho.isel(x=0), rho.isel(x=-1)
    at = lambda column, z: float(column.sel(z=z))
    assert at(left, -10e3) == AIR.params["rho"]            # sticky air
    assert at(left, -60e3) == WET_QUARTZ.params["rho"]     # crust
    assert at(left, -150e3) == DRY_OLIVINE.params["rho"]   # lithospheric mantle on the left
    assert at(right, -150e3) == WET_OLIVINE.params["rho"]  # asthenosphere below the thinner mantle on the right
    assert at(left, -250e3) == WET_OLIVINE.params["rho"]   # asthenosphere
    for key in tp.MATERIAL_PARAMETERS:
        assert key in builder.DTree.fields and not (builder.DTree.fields[key] == -999).any()


def test_export_interfaces_with_all_the_parameters(builder):
    builder.export_interfaces(export="none")
    lines = open(f"{builder.path}/interfaces.txt").read().splitlines()
    params = [line.split()[0] for line in lines[:len(tp.MATERIAL_PARAMETERS)]]
    assert params == ["C", "rho", "H", "A", "n", "Q", "V", "k", "weakening_seed", "cohesion_min", "cohesion_max",
                      "friction_angle_min", "friction_angle_max"]
    assert all(len(line.split()) == 1 + 4 for line in lines[:len(params)])  # one value per layer
    data = np.loadtxt(lines[len(params):])
    assert data.shape == (Nx, 3)  # one column per interface
    np.testing.assert_array_equal(data[:, 0], np.where(builder.x < Lx / 2, -160e3, -120e3))  # top of the asthenosphere
    np.testing.assert_array_equal(data[:, 2], -40e3)  # top of the crust


def test_params_keep_the_grid_of_the_builder(builder):
    builder.set_params(params_str="nx = 5\ntime_max = 1.0e6\nsolver = direct")
    assert builder.params["nx"] == Nx and builder.params["n_interfaces"] == 3
    assert builder.params["time_max"] == 1.0e6
    builder.export_params()
    text = open(f"{builder.path}/param.txt").read()
    assert "solver" in text and f"nx" in text and "dimensions" in text


def test_export_fields(builder):
    builder.create_temperature_field(_basic_temperature(builder))
    builder.export_field("temperature")
    T = np.loadtxt(f"{builder.path}/input_temperature_0.txt")
    assert T.shape == (Nx * Nz,)


# VelocityFieldBuilder ############################################################

def test_velocity_builder_conserves_the_flux(builder):
    with quiet():
        v = VelocityFieldBuilder(builder, units="cm/y")
        v.set_region("right", [-120e3, -40e3], vx=1.0)
        assert v.integrate_normal_components() > 0
        v.conservate_veloc("bot")
        assert abs(v.integrate_normal_components()) < 1e-12
        builder.create_velocity_field(v)
    vx = builder.DTree.fields.vx
    assert float(vx.sel(x=Lx, z=-80e3)) == pytest.approx(0.01 / tp.SEC_PER_YEAR)  # 1 cm/yr in m/s
    assert float(abs(vx.isel(x=slice(1, -1), z=slice(1, -1))).max()) == 0  # only the boundaries


@pytest.mark.parametrize("regime, moving_boundary", [("extension", "right"), ("extension", "left"),
                                                     ("extension", "symmetric"), ("compression", "right"),
                                                     ("compression", "left")])
def test_kinematic_regimes_conserve_the_flux(builder, regime, moving_boundary):
    with quiet():
        v = VelocityFieldBuilder(builder, units="cm/y")
        v.set_kinematic_regime(regime, 1.0, moving_boundary=moving_boundary, h_const=140e3)
        assert abs(v.integrate_normal_components()) < 1e-15
    left, right = v.velocs["left"].vx, v.velocs["right"].vx
    one = 0.01 / tp.SEC_PER_YEAR  # 1 cm/yr
    lithosphere = dict(z=-100e3)  # constant velocity (surface at -40 km, h_const = 140 km)
    air = dict(z=0.0)  # more than air_gap (10 km) above the surface: no velocity
    assert float(left.sel(**air)) == 0 and float(right.sel(**air)) == 0
    if moving_boundary == "symmetric":
        np.testing.assert_allclose(right.values, -left.values)
        assert float(right.sel(**lithosphere)) == pytest.approx(one / 2)
    elif regime == "extension":
        moving = right if moving_boundary == "right" else left
        sign = 1 if moving_boundary == "right" else -1
        assert float(moving.sel(**lithosphere)) == pytest.approx(sign * one)  # outwards
    else:
        moving, fixed = (right, left) if moving_boundary == "right" else (left, right)
        sign = -1 if moving_boundary == "right" else 1
        assert float(moving.sel(**lithosphere)) == pytest.approx(sign * one)  # inwards
        assert np.all(fixed.values == 0)


def test_symmetric_compression_is_an_error(builder):
    with quiet(), pytest.raises(ValueError, match="symmetric"):
        VelocityFieldBuilder(builder).set_kinematic_regime("compression", 1.0, moving_boundary="symmetric")


def test_linear_velocity(builder):
    with quiet():
        v = VelocityFieldBuilder(builder, units="m/s")
        v.linear_velocity("left", "vx", clims=[-Lz, 0], vlims=[1.0, 0.0])
    np.testing.assert_allclose(v.velocs["left"].vx.values, -builder.z / Lz, atol=1e-12)


# TemperatureFieldBuilder ##########################################################

ADIABAT = lambda depth, t_pot=1350.0: t_pot / np.exp(-10.0 * 3.28e-5 * depth / 1250.0)  # depth below the surface


def _basic_temperature(builder, t_pot=1350.0):
    """LAB: top of the asthenosphere (ID 0); surface: top of the crust (ID 2, the default)."""
    t = TemperatureFieldBuilder(builder, LAB_interface=0)
    t.apply_linear_temperature(t_pot=t_pot)
    return t


def test_temperature_builder_lithosphere(builder):
    t = TemperatureFieldBuilder(builder, LAB_interface=0)
    np.testing.assert_array_equal(t.surface_interface.values, -40e3)
    np.testing.assert_array_equal(t.lithosphere_thickness.values, np.where(builder.x < Lx / 2, 120e3, 80e3))


def test_linear_temperature(builder):
    T = _basic_temperature(builder).temperature
    assert float(T.min()) == 0.0
    assert float(T.sel(z=0).max()) == 0.0  # sticky air
    deep = T.sel(z=-250e3)  # asthenosphere: adiabat
    np.testing.assert_allclose(deep.values, ADIABAT(210e3))


def test_steady_columns_follow_the_lithosphere(builder):
    t = _basic_temperature(builder).solve_heat_1D()
    T = t.temperature
    np.testing.assert_allclose(T.sel(z=-40e3).values, 0.0, atol=1e-9)  # surface
    assert float(T.sel(z=-100e3).isel(x=-1)) > float(T.sel(z=-100e3).isel(x=0))  # thinner lithosphere is hotter
    builder.create_temperature_field(t)
    assert builder.DTree.fields.temperature.attrs["t_potential"] == 1350.0


def test_steady_2D_conducts_sideways(builder):
    T1 = _basic_temperature(builder).solve_heat_1D().temperature
    T2 = _basic_temperature(builder)
    T2.solve_heat_2D()
    T2 = T2.temperature
    np.testing.assert_allclose(T2.sel(z=-250e3).values, ADIABAT(210e3))  # asthenosphere kept
    np.testing.assert_allclose(T2.sel(z=0).values, 0.0)  # sticky air kept
    # far from the step of the LAB both solutions agree; at the step the cold side is heated
    far = T2.sel(z=-100e3).isel(x=0)
    assert float(far) == pytest.approx(float(T1.sel(z=-100e3).isel(x=0)), abs=1.0)
    near = dict(z=-100e3, x=Lx / 2 - 2e4)
    assert float(T2.sel(**near)) > float(T1.sel(**near)) + 10.0


def test_solve_heat_column_propagates_one_column(builder):
    t = _basic_temperature(builder)
    t.solve_heat_column(x_inx=0)
    np.testing.assert_array_equal(t.temperature.isel(x=0).values, t.temperature.isel(x=-1).values)


def test_yield_strength_envelope(builder):
    builder.create_temperature_field(_basic_temperature(builder))
    builder.calculate_YSE(sr=1e-15)
    fields = builder.DTree.fields
    for name in ["litho_P", "plastic_YSE", "viscous_YSE", "YSE", "YSE_soft"]:
        assert name in fields
    assert (fields.litho_P >= 0).all()
    assert (fields.YSE <= fields.plastic_YSE).all() and (fields.YSE <= fields.viscous_YSE).all()
    assert (fields.YSE_soft <= fields.YSE + 1e-6).all()  # softening lowers the strength
