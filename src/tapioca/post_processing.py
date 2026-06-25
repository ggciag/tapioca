import os, gc, json, sys, glob, pymp
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import RegularGridInterpolator

import numpy as np
import pandas as pd
import xarray as xr
from xarray import DataTree, register_datatree_accessor

from tapioca.scenClasses import MandyocScen

from ._variables import VARS_TYPES, SEC_PER_YEAR, VARIABLES_LIST

__all__ = ["scenPostProcessing"]

# Post processing objects and functions
@register_datatree_accessor("postproc")
class scenPostProcessing:
    """
    DataTree accessor for calculating common post-processing tasks and routines.

    This class registers under the ``postproc`` namespace for ``MandyocScen`` class,
    providing methods to compute physical properties directly on the 
    DataTree structure.

    Parameters
    ----------
    mandyocScen_DTree : xarray.DataTree
        The loaded scenario data tree containing the mesh and variables.
    varsTypes : dict, optional
        A dictionary defining the expected data types for calculated variables.
        Default is the global ``VARS_TYPES``.
    """

    def __init__(self, mandyocScen_DTree:xr.Datatree, varsTypes=VARS_TYPES):
        
        self.scenario = mandyocScen_DTree
        self.varsTypes = varsTypes
        
        return None
    
    def DeviatoricStressTensor(self, components='all', J2=True, mesh_upscaled=False,
                    units='MPa', export_nc=False):
        r"""
        Calculates the deviatoric stress tensor and second invariant (J2).

        This method computes the components of the deviatoric stress tensor 
        optimized for incompressible fluids (:math:`\nabla \cdot \mathbf{v} = 0`). 
        The calculated variables are appended directly to the ``'/mesh/original'`` 
        node of the DataTree.

        Parameters
        ----------
        components : list of str; or str, optional
            The specific tensor components to calculate. Currently defaults to 'all'.
            (Future improvement: filtering specific components for the 3D version).
        J2 : bool, optional
            If True, calculates the second invariant of the deviatoric stress 
            tensor (J2) and appends it to the mesh. 
            
            Default is True.
        mesh_upscaled : bool, optional
            If True, calculates the tensor using the upscaled lithology field 
            rather than the original mesh. (Future improvement). 
            
            Default is False.
        units : str, optional
            The physical units for the calculated stress. Supports 'MPa' 
            (scales the output by $10^6$) or standard Pascals. 
            
            Default is 'MPa'.
        export_nc : bool, optional
            If True, exports the resulting deviatoric stress tensor as a 
            NetCDF4 file. (Future improvement). 
            
            Default is False.

        
        Returns
        -------
        xarray.DataTree
            The updated scenario DataTree containing the new stress variables 
            (`tau_xx`, `tau_zz`, `tau_xz`, and optionally `tau_J2`).

        
        Notes
        -----
        **Signal Convention:**
        * Positive values indicate extension.
        * Negative values indicate compression.

        
        The deviatoric stress tensor :math:`\tau_{ij}` (or :math:`\sigma_{ij}'`) is calculated 
        through the relation presented in Gerya(2019):

        .. math::
            \tau_{ij} = 2 \eta \dot{\epsilon}_{ij}
        

        where :math:`\eta` is the effective viscosity and :math:`\dot{\epsilon}_{ij}` is the strain rate 
        tensor. The strain rate is derived from the velocity field :math:`v`:

        .. math::
            \dot{\epsilon}_{ij} = \frac{1}{2} \left( \frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i} \right)
        
        """

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
