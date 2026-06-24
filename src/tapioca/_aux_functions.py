import os, gc, json, sys, glob, pymp

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

from pathlib import Path

from ._variables import SEC_PER_YEAR, VARIABLES_LIST

__all__ = ["read_params","read_data"]

# Old functions, mostly made to handle the old format of data management
# Could be useful for people using older versions of Mandyoc

def read_params(path_param:str)->dict:
    '''
    Loads the param.txt in the form of a dictionary. 
    Keys are the parameter names and the values are the parameter values.
    
    Parameters
    ----------
    path_param : str or Path
        Path to the param.txt file.
    
    Returns
    -------
    dict
        Dictionary containing the parameter (key) and its value (value).
    '''

    params_form = {}
    ptemp = ''
    with open(path_param,'r') as param:
        for line in param:
            line = line.strip()
            if len(line)==0:
                continue
            elif line[0] == "#":
                continue
            
            line = line.split('#')[0]
            line = line.replace(' ','')
            pv = line.split('=')
            params_form[pv[0].lower()] = pv[1]
            ptemp = ptemp + line+'\n'

    return params_form

def read_data(file: str, Nx: int, Nz: int, veloc:bool=False, surface:bool=False) -> np.array:
    '''
    Reads the data from a .txt file for a variable and returns it as a numpy array. 
    The function can also return the velocity components or surface data separately.

    Parameters
    ----------
    file : str or Path
        File to read (.txt)
    
    Nx : int
        Number of elements in the x direction
        
    Nz : int
        Number of elements in the z direction

    veloc : bool, optional
        If True, returns the velocity components separately.
        Default is False.
    
    surface : bool, optional
        If True, returns the surface data.
        Default is False.

    Returns
    -------
    np.array
        Variable loaded.

    Notes
    ----
    - If veloc is True, returns the velocity components separately.
    - If surface is True, returns the surface data.
    
    Otherwise, returns the data as a numpy array.
    '''

    #file = f'{"_".join(file.split("_")[:-1])}/{file}.txt'
    data = pd.read_csv(file, header=None, 
                       skiprows=2, comment='P')
    
    data = data.to_numpy()

    if not(veloc) and not(surface):
        data[np.abs(data) < 1.0e-200] = 0 #converter numeros pequenos e grandes
        data = np.reshape(data, (Nx,Nz), order='F') #(nx*nz,1) -> (nx,nz)
        data = data.T
        
    elif surface == True:
        return data
        
    else:
        vx = np.reshape(data[0::2], (Nx,Nz), order='F')
        vy = np.reshape(data[1::2], (Nx,Nz), order='F')
        data = (vx.T, vy.T)
    return data


#====== OLD/DEPRECATED FUNCTIONS ======

def _old_get_rank(cdir:str) -> int:
    '''Get the number of processes (ranks) used to run a
    scenario based on the steps directory.

    Notes
    -----
    This is an old function, kept to test somethings.

    Parameters
    ----------
    cdir : str
        Path to the scenario directory.

    Returns
    -------
    int
        Number of ranks used to run the scenario.

    '''

    return int(len(list(Path(f'{cdir}/steps').glob("step_0_*"))))

def _old_get_lasttime(cdir:str) -> int:
    '''
    Get the last time step of a scenario based on the time directory.

    Notes
    -----
    This is an old function, kept for compatibility.

    Parameters
    ----------
    cdir : str
        Path to the scenario directory.

    Returns
    -------
    int
        Last time step of the scenario.
    '''
    files = glob.glob(os.path.join(cdir,'time/time*'))
    times = []
    for f in files:
        times.append(int(f[:-4].split('_')[-1]))
    
    return np.max(times)