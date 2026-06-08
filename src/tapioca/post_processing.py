import os, gc, json, sys, glob, pymp
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import RegularGridInterpolator

import numpy as np
import pandas as pd
import xarray as xr
from xarray import DataTree, register_datatree_accessor

from ._variables import VARS_TYPES, SEC_PER_YEAR, VARIABLES_LIST

# Post processing objects and functions
@register_datatree_accessor("postproc")
class scenPostProcessing:

    def __init__(self, mandyocScen_DTree, varsTypes=VARS_TYPES):
        '''
        It is a class to calculate common post-processing variables on a loaded scenario
        '''
        
        self.scenario = mandyocScen_DTree
        self.varsTypes = varsTypes
        
        return None
    
    def DeviatoricStressTensor(self, components='all', J2=True, mesh_upscaled=False,
                    units='MPa', export_nc=False):
        '''
        Calculates the deviatoric stress tensor and second invariant (J2).
        Optimized for incompressible fluids (div(v) = 0).
        
        Signal convention:
           Positive: extension
           Negative: compression


        ---
        
        This function calculate the deviatoric stress tensor (\\tau or \\sigma') through the relation:
           $ \\tau_{ij} = 2 d \\epsilon_{ij} \\dt * \\nu $
        where ij are the coodinates components, d \\epsilon \\dt is the strain rate, and \nu is the viscosity.
        
        The strain rate is calculate using the velocity field (v_{ij}) by:
            d \\epsilon_{ij} \\dt = 1/2 (dv_{i}/dx_{j} + dv_{j}/dx_{i}) 
        
        ---

        components -> future improvement: to calculate only the desired components (useful for the 3D version)
        mesh_upscaled -> future improvement: to use a upscaled field rather than the original mesh
        export_nc -> future improvement: to export the deviatoric stress tensor as a NETCDF4 file
        '''

        dtype = self.varsTypes['deviatoric_stress']
        tau_factor = 1.0
        if units == 'MPa':
            tau_factor = 1e6

        
        scen = self.scenario
        mesh = scen.mesh['original'].to_dataset()

        vx = mesh['vx']
        vz = mesh['vy'] # current export uses vy instead vz, it must be rewrite for the 3D version
        visc = mesh['viscosity']

        # calculate the strain tensor components
        dvx_dx = vx.differentiate('x',edge_order=2)
        dvx_dz = vx.differentiate('z',edge_order=2)   #\\tau_{xz} = \\tau_{zx} 
        dvz_dx = vz.differentiate('x', edge_order=2)

        # calculate the deviatoric stress tensor components 
        self.scenario['/mesh/original']['tau_xx'] = (2 * visc * dvx_dx / tau_factor).astype(dtype)
        self.scenario['/mesh/original']['tau_zz'] = -self.scenario['/mesh/original']['tau_xx'] # considering fluid incompressibility: u_{i,i} = 0
        
        self.scenario['/mesh/original']['tau_xz'] = (visc * (dvx_dz + dvz_dx) / tau_factor).astype(dtype) 

        for comp in ['tau_xx','tau_zz','tau_xz']:
            self.scenario['/mesh/original'][comp].attrs['long_name'] = comp
            self.scenario['/mesh/original'][comp].attrs['units'] = units
            self.scenario['/mesh/original'][comp].attrs['description'] = f'component {comp} of the deviatoric stress tensor'
        
        self.scenario['/mesh/original']['tau_xz'].attrs['note'] = 'tau_xz = tau_zx'
        
        if J2==True:
            self.scenario['/mesh/original']['tau_J2'] = 1/2 * (self.scenario['/mesh/original']['tau_xx']**2 + self.scenario['/mesh/original']['tau_zz']**2) + self.scenario['/mesh/original']['tau_xz']**2
            self.scenario['/mesh/original']['tau_J2'].attrs['units'] = units
            self.scenario['/mesh/original']['tau_J2'].attrs['long_name'] = 'second invariant of the deviatoric stress tensor'
        
        if export_nc:
            pass # Implement xarray to_netcdf here

        
        return self.scenario
