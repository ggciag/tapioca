import numpy as np
from .generic_input_classes import MandyocBuilder

__all__=['get_interface']

def get_interface(scenario:MandyocBuilder, reference:str|int):
    """
    Get a scenario builder (`MandyocBuilder`) and return an interface 
    DataArray by id (int) or name (str).
    
    
    Parameters
    ----------
    scenario: MandyocBuilder


    Return
    ------
    xr.DataArray
    
    """
    ds = scenario.DTree.interfaces.ds

    if isinstance(reference, str):
        if reference not in ds.data_vars:
            raise KeyError(f"Interface '{reference}' not found.")
        return ds[reference]

    elif isinstance(reference, (int, np.integer)):
        for var in ds.data_vars:
            if int(ds[var].attrs.get('id', -999)) == int(reference):
                return ds[var]

