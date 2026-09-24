import os, gc, json, sys, glob

import numpy as np
import xarray as xr
import pandas as pd

from pathlib import Path

from ._variables import SEC_PER_YEAR, VARIABLES_LIST, INTERFACES_PARAMETERS

__all__ = ["read_params","read_data","ensure_directory_exists","export_compressed_dataset","export_compressed_datatree",
           "_numba_diffusion_loop",
           "_Druker_Prager_YS","_Byerlee_law_YS","_visc_dislocation_creep"]

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
    -----
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

def ensure_directory_exists(folder_path: str|Path):
    """
    Checks if a directory path exists and creates it if it does not. 
    Returns True if succeed. 

    Parameters
    ----------
    folder_path : str or pathlib.Path
        The path of the directory to check and create.

    """
    path = Path(folder_path)
    path.mkdir(parents=True, exist_ok=True)

    return True

def export_compressed_dataset(ds, filename, complevel=5):
    """
    Exports an xarray.Dataset to a compressed NetCDF4 file.
    
    Parameters
    ----------
    dt : xarray.Dataset
        The Dataset containing the groups.
    filename : str
        The output file path.
    complevel : int
        Compression level from 1 (fastest) to 9 (smallest file size). 
        Default is 5.
    """

    for var_name, var_data in ds.data_vars.items():
        var_data.encoding.update({
            'zlib': True, 
            'complevel': complevel
        })
            
    # Export the dataset to a NetCDF4 file
    ds.to_netcdf(filename, engine="h5netcdf")

    print(f"Successfully exported compressed Dataset to {filename}")

def export_compressed_datatree(dt, filename, complevel=5):
    """
    Exports an xarray.DataTree to a compressed NetCDF4 file.
    
    Parameters
    ----------
    dt : xarray.DataTree
        The DataTree containing the groups.
    filename : str
        The output file path.
    complevel : int
        Compression level from 1 (fastest) to 9 (smallest file size). 
        Default is 5.
    """
    
    # Iterate through every group/node in the DataTree
    for node in dt.subtree:
        # Iterate every data variable in the current group
        for var_name, var_data in node.data_vars.items():
            var_data.encoding.update({
                'zlib': True, 
                'complevel': complevel
            })
            
    # Export the DataTree to a single NetCDF4 file
    dt.to_netcdf(filename, engine="h5netcdf")


from numba import njit
@njit(fastmath=True)
def _numba_diffusion_loop(T, kappa, H, c_cap, dx, dz, dt, num_steps, cond):
    """JIT-compiled loop for performance gains
    It is a function to iterate the Euler Foward Approach in the heat diffusion equation.
    """
    for step in range(num_steps):
        T_old = np.copy(T)
        T_new = np.copy(T)
        
        dT_dx = (T[2:, 1:-1] - T[:-2, 1:-1]) / (2.0 * dx)
        dT_dz = (T[1:-1, 2:] - T[1:-1, :-2]) / (2.0 * dz)
        
        dK_dx = (kappa[2:, 1:-1] - kappa[:-2, 1:-1]) / (2.0 * dx)
        dK_dz = (kappa[1:-1, 2:] - kappa[1:-1, :-2]) / (2.0 * dz)
        
        d2T_dx2 = (T[2:, 1:-1] - 2.0 * T[1:-1, 1:-1] + T[:-2, 1:-1]) / (dx**2)
        d2T_dz2 = (T[1:-1, 2:] - 2.0 * T[1:-1, 1:-1] + T[1:-1, :-2]) / (dz**2)
        
        diffusion_x = kappa[1:-1, 1:-1] * d2T_dx2 + dK_dx * dT_dx
        diffusion_z = kappa[1:-1, 1:-1] * d2T_dz2 + dK_dz * dT_dz
        
        T_new[1:-1, 1:-1] = T[1:-1, 1:-1] + dt * (diffusion_x + diffusion_z + H[1:-1, 1:-1] / c_cap)
        
        # Boundaries
        T_new[0, :] = T_new[1, :]    
        T_new[-1, :] = T_new[-2, :]  
        T_new[:, 0] = T[:, 0]
        T_new[:, -1] = T[:, -1]
        
        # Apply mask
        T = np.where(cond, T_new, T)
            
    return T

def _Druker_Prager_YS(c:float, Phi:float, P:float):
    '''
    Gives the plastic Yield Stress (in Pa) according to the Druker-Prager criterion:

    Tau_yield = c * cos(phi) + P * sin(phi)

    where c if the internal cohesion (in Pa), phi is the internal angle of friction (in degrees), 
    and P is the pressure (in Pa).
    '''
    rad = np.pi/180
    return c * np.cos(Phi*rad) + P * np.sin(Phi*rad)

def _Byerlee_law_YS(c:float, mu:float, P:float):
    '''
    Gives the plastic Yield Stress (in Pa) according to the Byerlee Law:

    Tau_yield = c + mu * P

    where c if the internal cohesion (in Pa), mu is the friction coefficient (dimensionless), 
    and P is the pressure (in Pa).
    '''
    return c + mu * P

def _visc_dislocation_creep(strain_rate:float, A:float, C:float, n:float, Q:float, V:float, P:float, R:float, T:float):
    '''
    Return the viscosity of a deslocation creep flow for a given strain rate. 
    
    It is calculated using an Arrhenius-type constitutive equation implemented in mandyoc:
        
        https://ggciag.github.io/mandyoc/files/implementation.html#rheology
    '''

    return C * A**(-1./n) * strain_rate**((1.0-n)/n)*np.exp((Q + V*P)/(n*R*T))

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