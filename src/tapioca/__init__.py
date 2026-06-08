import warnings
import xarray as xr
import numpy as np

__version__ = "0.1.0"

# 1. Import variables
from ._variables import *

# 2. Import Scenario Classes
from .scenClasses import *

with warnings.catch_warnings(): #filtering warning for the preexisting attribute of the datatree accessor
    warnings.filterwarnings("ignore", message=".*name 'postproc'.*")
    # 3. Import post-processing classes and functions
    from .post_processing import *

# 5. Import plotting functions


# 4. Import utils functions
from ._aux_functions import *

