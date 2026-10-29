"""Tests of the checks of the input parameters (tapioca.inputs.checks)."""

import pytest

from tapioca.inputs.checks import InputChecks, InputError, SVGChecks, SVGInputError

from conftest import build_scenario


@pytest.fixture
def checks():
    return InputChecks()


def test_choice(checks):
    assert checks._check_choice("solver", "2d", ["2D", "columns"]) == "2D"
    assert checks._check_choice("solver", "1D", ["2D", "columns"], aliases={"1D": "columns"}) == "columns"
    with pytest.raises(InputError, match="Unknown solver 'x'. Use 2D, columns."):
        checks._check_choice("solver", "x", ["2D", "columns"])


def test_number(checks):
    assert checks._check_number("v", "1.5") == 1.5
    assert checks._check_number("t", 0, 0) == 0
    with pytest.raises(InputError, match="greater than 0"):
        checks._check_number("t", 0, 0, strict=True)
    with pytest.raises(InputError, match="must be a number"):
        checks._check_number("v", "fast")


def test_path(checks, tmp_path):
    (tmp_path / "a.txt").write_text("")
    assert checks._check_path("file", "a.txt", base=str(tmp_path)) == str(tmp_path / "a.txt")
    with pytest.raises(InputError, match="file not found"):
        checks._check_path("file", "b.txt", base=str(tmp_path))


def test_keys_and_domain(checks):
    checks._check_keys("setting", ["a"], ["a", "b"])
    with pytest.raises(InputError, match="Unknown setting: c"):
        checks._check_keys("setting", ["a", "c"], ["a", "b"])
    with pytest.raises(InputError, match="Nx must be at least 2"):
        checks._check_domain(1e3, 1e3, 1, 10)


def test_svg_checks_raise_svg_errors():
    assert issubclass(SVGInputError, InputError) and issubclass(SVGInputError, ValueError)
    with pytest.raises(SVGInputError):
        SVGChecks()._check_on_off("magmatism", "yes")


def test_unknown_setting_of_the_scenario():
    with pytest.raises(SVGInputError, match="Unknown setting: regim"):
        build_scenario("rifting", regim="compression")


@pytest.mark.parametrize("settings, message", [({"quiescence_time": -1}, "quiescence_time must be at least 0"),
                                               ({"magmatism": "yes"}, "Unknown magmatism"),
                                               ({"regime": "drift"}, "Unknown regime")])
def test_invalid_settings(settings, message):
    with pytest.raises(SVGInputError, match=message):
        build_scenario("rifting", **settings)


def test_settings_with_units():
    scenario = build_scenario("accordion", quiescence_time="300 Myr", velocity="5 mm/yr")
    assert scenario.quiescence_time == 300.0 and scenario.velocity == 0.5
