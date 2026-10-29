"""
Input builders for Mandyoc scenarios drawn in SVG files (e.g. made in Inkscape).

- materials: Mandyoc layers (MandyocLayer) and their rheologies;
- reader: geometry and materials of a scenario drawn in an SVG file;
- scenario: SVGScenario, which builds all the input files of a scenario from an SVG
  (interfaces, parameters, initial temperature and velocity for rifting, compression,
  accordion and quiescence experiments).
"""

# getting the same version of the whole module
from ...__init__ import __version__

from .reader import *
from .materials import *
from .scenario import *
