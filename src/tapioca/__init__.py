"""
Tapioca-Mandyoc: Input and Output management and post-processing for Mandyoc geodynamic models.
"""

from warnings import catch_warnings, filterwarnings

__version__ = "0.1.0"

# 1. Import variables
from ._variables import SEC_PER_YEAR, VARIABLES_LIST, VARS_TYPES

# 2. Import Scenario Classes
from .scenClasses import MandyocScen

# 5. Import plotting functions


# 4. Import utils functions
from ._aux_functions import *

