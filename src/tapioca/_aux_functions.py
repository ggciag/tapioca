import os, gc, json, sys, glob

import numpy as np
import xarray as xr
import pandas as pd

import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numba import jit

from pathlib import Path

from ._variables import SEC_PER_YEAR, VARIABLES_LIST, INTERFACES_PARAMETERS,PARAMS_TYPES

__all__ = ["read_params","evaluate_params","read_data","ensure_directory_exists","export_compressed_dataset","export_compressed_datatree",
           "_numba_diffusion_loop","build_diagonals_2D","spsolve_heat_col",
           "_Druker_Prager_YS","_Byerlee_law_YS","_visc_dislocation_creep"]

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
            params_form[pv[0]] = PARAMS_TYPES[pv[0]](pv[1])
            ptemp = ptemp + line+'\n'

    return params_form

def evaluate_params(param_file:str='', params_str:str='', **params_dict)->dict:
    temp_param = {}

    if len(param_file) > 1:
        temp_param = read_params(param_file)
        
    elif len(params_str) > 1:
        for line in params_str.split('\n'):
            line = line.strip()

            if len(line)==0:
                continue
            elif line[0] == "#":
                continue

            line = line.split('#')[0]
            line = line.replace(' ','')
            pv = line.split('=')
            temp_param[pv[0]] = PARAMS_TYPES[pv[0]](pv[1])
    else:
        for k,v in params_dict.items():
            temp_param.update({k:PARAMS_TYPES[k](v)})

    return temp_param

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


@jit(nopython=True,fastmath=True)
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

def spsolve_heat_col(T, k_nodes, H, Tsurf, Tlab, dz):
    """
    This function solves the 1D steady state heat diffusion equation on a single vertical 
    column. The formulation is:

    dT/dt = d/dz (k dT/dz) + H 

    Parameters
    ----------
    T : ndarray
        1D array of the temperature profile. In this steady-state solver, its values 
        are not used for computation, just dictates the number of nodes (N = len(T)).
    k_nodes : ndarray
        1D array of thermal conductivities at the nodes. 
        Must be of length N.
    H : ndarray
        1D array of volumetric radiogenic heat production (W/m^3) at each node. 
        Must be of length N.
    Tsurf : float
        Temperature boundary condition. Note: In the current implementation, this 
        is applied to the bottom node (index -1).
    Tlab : float
        Temperature boundary condition. Note: In the current implementation, this 
        is applied to the top node (index 0).
    dz : float
        Vertical spacing between nodes (meters).

    Returns
    -------
    ndarray
        1D array of length N containing the solved steady-state temperature profile.
    """

    k = 2* (k_nodes[:-1] * k_nodes[1:]) / (k_nodes[1:] + k_nodes[:-1])
    coeffs = k / (dz **2)
    N = len(T)
    
    # k_{i+1/2} / dz**2
    lower_diag = np.copy(coeffs) 
    upper_diag = np.copy(coeffs)
    
    main_diag = np.zeros(N)
    d = -np.copy(H) 

    # -( (k_{i+1/2} + k_{i-1/2})/dz**2 )
    main_diag[1:-1] = -(coeffs[:-1]+coeffs[1:])

    # Boundary conditions
    main_diag[0] = 1.0
    upper_diag[0] = 0.0
    d[0] = Tlab # T surface

    main_diag[-1] = 1.0
    lower_diag[-1] = 0.0
    d[-1] = Tsurf

    # Assemble the sparse matrix (-1: lower diagonal, 0:main diag, 1:upper diagonal)
    A = sp.diags([lower_diag, main_diag, upper_diag], [-1, 0, 1], format='csc')
    T = spla.spsolve(A, d) # A*T=d
    return T

@jit(nopython=True,cache=True)
def build_diagonals_2D(Nz, Nx, dz, dx, k_2D, H_2D, T_top,ZZ,Tmesh,lab_surf,bot_bc='adiabatic'):
    """
    Fast matrix assembly using Numba. 
    This functions is made to solve 2D steady state heat equation.
    i=0 is bottom, i=Nz-1 is top. j=0 is left, j=Nx-1 is right.
    """
    N = Nz * Nx
    
    # Pre-allocate the 5 diagonals
    diag_main = np.zeros(N)
    diag_left = np.zeros(N - 1)
    diag_right = np.zeros(N - 1)
    diag_down = np.zeros(N - Nx)
    diag_up = np.zeros(N - Nx)
    
    # Pre-allocate Right-Hand Side (b)
    b = np.zeros(N)
    
    for i in range(Nz):
        for j in range(Nx):
            idx = i * Nx + j

            # --- 1. Apply Top and Bottom Boundary Conditions ---
            # check air
            if ZZ[i,j] >= 0:
                diag_main[idx] = 1.0
                b[idx] = T_top
                continue

            if bot_bc.lower() == 'adiabatic': # keeps pre set temperature
                # check asthenosphere
                if ZZ[i,j] <= -lab_surf[j]:
                    diag_main[idx] = 1.0
                    b[idx] = Tmesh[i,j]
                    continue
            
            else:
                if i == 0: # fix only the bottom of the model
                    diag_main[idx] = 1.0
                    b[idx] = Tmesh[i,j]
                    continue

            # air (i = Nz-1), in case there is no sticky air layer
            if i == Nz - 1:
                diag_main[idx] = 1.0
                b[idx] = T_top                 # fallback: treat top row as surface
                continue
            # --- 2. Interior and Side Boundaries ---
            
            # Vertical harmonic means (Always exist for interior nodes)
            k_dn = 2.0 * k_2D[i, j] * k_2D[i-1, j] / (k_2D[i, j] + k_2D[i-1, j] + 1e-15)
            k_up = 2.0 * k_2D[i, j] * k_2D[i+1, j] / (k_2D[i, j] + k_2D[i+1, j] + 1e-15)
            
            c_dn = k_dn / (dz**2)
            c_up = k_up / (dz**2)
            c_lt = 0.0
            c_rt = 0.0
            
            # Horizontal harmonic means (Respecting Neumann side boundaries)
            if j > 0: # Has left neighbor
                k_lt = 2.0 * k_2D[i, j] * k_2D[i, j-1] / (k_2D[i, j] + k_2D[i, j-1] + 1e-15)
                c_lt = k_lt / (dx**2)
                
            if j < Nx - 1: # Has right neighbor
                k_rt = 2.0 * k_2D[i, j] * k_2D[i, j+1] / (k_2D[i, j] + k_2D[i, j+1] + 1e-15)
                c_rt = k_rt / (dx**2)
            
            # --- Normalize the equation by the main diagonal ---
            norm = c_dn + c_up + c_lt + c_rt
            
            if j > 0:
                diag_left[idx - 1] = c_lt / norm
            if j < Nx - 1:
                diag_right[idx] = c_rt / norm
                
            # Place vertical coefficients scaled
            diag_down[idx - Nx] = c_dn / norm
            diag_up[idx] = c_up / norm
            
            # Place main diagonal and RHS scaled
            diag_main[idx] = -1.0
            b[idx] = -H_2D[i, j] / norm
            
    return diag_down, diag_left, diag_main, diag_right, diag_up, b

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


def get_interface(interfaces_ds, reference:str|int):
    """
    Get a scenario builder (`MandyocBuilder`) and return an interface 
    DataArray by id (int) or name (str).
    
    
    Parameters
    ----------
    interfaces_ds: xr.Dataset
        The dataset with all interfaces

    reference: str or int
        The reference to get the interface. If is a string, gets by name, 
        if is a int, gets by ID.
    Return
    ------
    xr.DataArray
        The selected interface (Nx).
    """

    if isinstance(reference, str):
        if reference not in interfaces_ds.data_vars:
            raise KeyError(f"Interface '{reference}' not found.")
        return interfaces_ds[reference]

    elif isinstance(reference, (int, np.integer)):
        for var in interfaces_ds.data_vars:
            if int(interfaces_ds[var].attrs.get('id', -999)) == int(reference):
                return interfaces_ds[var]
            
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