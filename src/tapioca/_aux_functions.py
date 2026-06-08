import os, gc, json, sys, glob, pymp

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

from pathlib import Path

from ._variables import *







# Old functions, mostly made to handle the old format of data management
# Could be useful for people using older versions of Mandyoc

def read_params(path_param):
    '''
    [Old function]
    This functions loads param.txt in the form of a dictionary. 
    Keys are the parameter names and the values are the parameter values. 
    
    args:
        path_param: path to the param.txt file of a scenario.
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

def get_rank(cdir):
    '''
    [Old function]
    This function returns the number of processes (ranks) used to run a scenario based on the steps directory.

    args:
        cdir: path to the scenario directory.

    return:
        int: number of ranks used to run the scenario.
    '''

    return int(len(list(Path(f'{cdir}/steps').glob("step_0_*"))))

def get_lasttime(cdir):
    '''
    [Old function]
    This function returns the last time step of a scenario based on the time directory.
    
    args:
        cdir: path to the scenario directory.

    return:
        int: last time step of the scenario.
    '''
    files = glob.glob(os.path.join(cdir,'time/time*'))
    times = []
    for f in files:
        times.append(int(f[:-4].split('_')[-1]))
    
    return np.max(times)

def read_data(file, Nx, Nz, veloc=False, surface=False):
    '''
    [Old function]
    This function reads the data from .txt files for a variable and returns it as a numpy array. 
    The function can also return the velocity components or surface data separately.

    args:
        file: .txt file to read.
        veloc: boolean, if True, returns the velocity components separately.
        surface: boolean, if True, returns the surface data.

        Nx: number of grid points in the x direction. [can be get from the params]
        Nz: number of grid points in the z direction. [can be get from the params]

    return:
        numpy array: variable loaded.
    '''

    file = f'{"_".join(file.split("_")[:-1])}/{file}.txt'
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