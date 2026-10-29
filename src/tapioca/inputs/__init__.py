
"""
- MandyocBuilder, VelocityFieldBuilder and TemperatureFieldBuilder: generic builders of the
  interfaces, material fields, velocity and temperature of a scenario, used from scripts;
- checks: InputChecks (and SVGChecks), standardized checks of the input parameters;

- materials: Mandyoc layers (MandyocLayer) and their rheologies;
- svg: geometry and materials of a scenario drawn in an SVG file (e.g. made in Inkscape);
- scenario: SVGScenario, which builds all the input files of a scenario from an SVG
  (interfaces, parameters, initial temperature and velocity for rifting, compression,
  accordion and quiescence experiments).
"""

from ..__init__ import __version__

# getting the generic input classes

from .generic_input_classes import * 
from .checks import InputError, SVGInputError

from . import svg

__all__ = ["MandyocBuilder", "VelocityFieldBuilder", "TemperatureFieldBuilder", "InputError", "SVGInputError", "svg"]


