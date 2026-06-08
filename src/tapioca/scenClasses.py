import os, gc, json, sys, glob, pymp
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import RegularGridInterpolator

import numpy as np
import pandas as pd
import xarray as xr
from xarray import DataTree, register_datatree_accessor

from ._variables import VARS_TYPES, VARIABLES_LIST

#Mandyoc Scenario class
class MandyocScen:
    
    def __init__(self, path, variables=['density'], name=None,
                 load_lithology=False, load_surface=False, load_particles=False,
                 particles_file='particles_trajectories.nc',
                 xlimits=None, zlimits=None, tlimits=None, # ylimits should be implemented for the 3D version
                 thick_air=40e3,
                 chunks_vars={"x": 'auto', "z": 'auto', 'time': "auto"},
                 filter_air=True, #only relevant if load particles
                 verbose=False):

        # Setting directories and scen name
        self.path = Path(path)
        self.verbose = verbose
        if isinstance(name, str): self.name = name
        else: self.name = self.path.name

        if self.verbose:
            print(f'Scenario at: {self.path}')
            print(f'Scenario name: {self.name}')
        
        # Handle Variables and Extract Metadata
        if isinstance(variables, str): variables = [variables]
        elif not isinstance(variables, (list, tuple, np.ndarray)): variables = list(variables)
        
        if len(variables) == 0:
            raise ValueError("The 'variables' list cannot be empty. Need at least one variable to extract metadata.")
        
        self.get_scenarioData(variables[0])
        
        self.xlimits = xlimits if xlimits is not None else [self.XMIN, self.XMAX]
        self.zlimits = zlimits if zlimits is not None else [self.ZMIN, self.ZMAX]
        self.tlimits = tlimits if tlimits is not None else [self.TMIN, self.TMAX]
        self.thick_air = thick_air #m

        if self.verbose:
            print(f'x limits: {self.xlimits}')
            print(f'z limits: {self.zlimits}')
            print(f'time limits: {self.tlimits}')
        
        self.z_corrected = False
        self.particles_loaded = False
        
        # Initialize the empty DataTree
        self.DTree = DataTree()

        # Passing some metadata to the DataTree
        self.DTree.attrs['name'] = self.name
        self.DTree.attrs['xlimits'] = self.xlimits
        self.DTree.attrs['zlimits'] = self.zlimits
        self.DTree.attrs['tlimits'] = self.tlimits
        
        # Reading and storing standard variables (Nx,Nz)
        standard_datasets = []
        for var in variables:
            ds = self._load_spatial_var(var, chunks=chunks_vars)
            if ds is not None:
                standard_datasets.append(ds)
        
        if standard_datasets:
            self.DTree['/mesh/original'] = xr.merge(standard_datasets)

        if self.verbose: print(f"Variables loaded: {' '.join(variables)}")
        
        # Reading and storing the lithology (upscaled mesh)
        if load_lithology==True:
            ds_litho = self._load_spatial_var('lithology', chunks=chunks_vars)
            if ds_litho is not None:
                self.DTree['/mesh/upscaled'] = ds_litho

        if self.verbose: print("Lithology loaded")
            
        # Reading the surface/topography (only X dimension)
        if load_surface==True:
            
            ds_surf = self._load_spatial_var('surface', chunks=chunks_vars)
            if ds_surf is not None:
                self.DTree['/surface/topography'] = ds_surf

        if self.verbose: print("Surface loaded")
            
        # Loading particles
        if load_particles:
            self._load_particles(particles_file, chunks={'id': 'auto'}, filter_air=filter_air)
            
        if self.verbose: print(f"Particles [{particles_file}] loaded")
            
        #self.original_particles = None   #whole particles dataset (can be replaced for subsets)
        self.selected_particles = None   #current selected particles
        self.particles = {}  #dictionary with particle selections
        
        return None
    
    def get_scenarioData(self, var='density'):
        file_path = self.path / f"{var}.nc"
            
        if not file_path.exists():
            raise FileNotFoundError(f"Cannot extract metadata. File not found: {file_path}")

        # Open lazily and safely close automatically using 'with'
        with xr.open_dataset(file_path) as ds:
            # Getting dimensions, maximum and minimum from a base netcdf
            self.Nx = ds.sizes['x']
            self.Nz = ds.sizes['z']
                        
            self.XMAX = ds.x.max().item()
            self.XMIN = ds.x.min().item()
            self.ZMAX = ds.z.max().item()
            self.ZMIN = ds.z.min().item()
            self.TMAX = ds.time.max().item()
            self.TMIN = ds.time.min().item()

        return True
    
    def correctZcoord(self, factor=None):
        if self.z_corrected:
            if self.verbose: print("Z was already corrected")
            return False
        
        if factor is None:
            factor = self.ZMAX - self.thick_air
        
        # Traverse specific tree nodes and apply Z-correction if 'z' exists
        structs_to_check = ['/mesh', '/surface','/particles']
        
        mesh_nodes = list(self.DTree['/mesh'].children) # Mesh data (Eulerian grids)
        surface_nodes = list(self.DTree['/surface'].children) # Surface/string data (only X dimension)
        particles_nodes = list(self.DTree['/particles'].children) # Lagrangian particles 
        
        for node in mesh_nodes:
            node_ds = self.DTree['/mesh'][node].ds
            self.DTree['/mesh'][node] = node_ds.assign_coords(z=(node_ds['z'] - factor))
        
        # Add a If condition for when the used does not load all data types
        for node in surface_nodes:    #O script em julia está exportando com a superfície corrigida, mudar para exportar com o dado ORIGINAL
            node_ds = self.DTree['/surface'][node].ds
            self.DTree['/surface'][node] = node_ds.assign(surface=(node_ds['surface'] - factor))

        for node in particles_nodes:
            node_ds = self.DTree['/particles'][node].ds
            self.DTree['/particles'][node] = node_ds.assign(z=(node_ds['z'] - factor))

        self.zlimits = [self.zlimits[0] - factor, self.zlimits[1] - factor]
        self.ZMAX -= factor
        self.ZMIN -= factor
        
        if self.verbose:
            print(f"Z coordinate corrected by subtracting {factor} m")
            print(f"New z limits: {self.zlimits}")
            
        self.z_corrected = True
        return True

    def _load_spatial_var(self, variable, chunks={}):
        """
        Internal function to load Eulerian grids and slice them to xlimits,zlimits
        """
        
        file_path = self.path / f"{variable}.nc"
        
        if not file_path.exists():
            print(f"Warning: {file_path} not found. Skipping.")
            return None
            
        v = xr.open_dataset(file_path, chunks=chunks)
        
        if 'z' in v.coords:
            v = v.sortby('z')
            v = v.sel(
                x=slice(self.xlimits[0], self.xlimits[-1]),
                z=slice(self.zlimits[0], self.zlimits[-1]),
                time=slice(self.tlimits[0], self.tlimits[-1])
            )
        else:
            # Surface or 2D variables
            v = v.sel(x=slice(self.xlimits[0], self.xlimits[-1]),
                     time=slice(self.tlimits[0], self.tlimits[-1]))
            
        return v
    
    def _load_particles(self, name, chunks={'id': 'auto'}, filter_air=True, air_layer=None):
        """
        Internal function to load particles dataset 
        """        
        file_path = self.path / name
        
        if not file_path.exists():
            print(f"Warning: Particles file {file_path} not found. Skipping.")
            return False
            
        particles = xr.open_dataset(file_path, chunks=chunks)
        particles = particles.sel(time=slice(self.tlimits[0], self.tlimits[-1]))
        
        if filter_air and 'layer' in particles.data_vars:
            if isinstance(air_layer, int): air = air_layer
            else: air = int(particles.layer.max())
                
            cond = particles != air
            particles = particles.where(cond)
            if verbose: print(f'Air particles filtered [{air}]')
            
        self.DTree['/particles/original'] = particles
        self.particles_loaded = True
        return True
    
    
    def _apply_selection(self, valid_ids, replace_original=True, selected_name='',
                                    selection_name=''):
        """
        Internal function to select a dataset by the IDs
        """
        
        if replace_original==True:
            selection_name = 'original'
        elif selection_name == '':
            selection_name == 'selected'
        
        source_ds = self.DTree.particles[selected_name].ds
    
        self.DTree.particles[selection_name] = source_ds.sel(id=valid_ids)
    
        gc.collect()
        return None

    def _get_pts(self, select_original=True, selected_name=''):
        if select_original==True: selected_name = 'original'
        elif selected_name == '': selected_name = 'selected'
        pts = self.DTree.particles[selected_name].ds
        return pts, selected_name
    
    #Selecting particles
    def selectParticles_bytimerange(self, timerange, select_original=True, selected_name='',
                                    selection_name='', replace_original=False):
        '''
        select particles that appeared (e.g. sedimented) within the specified time range
        timerange : array-like = [tmin, tmax]
        '''
        # apply the support selection function
        pts, selected_name = self._get_pts(select_original, selected_name)
        
        # get IDs of the time "i" and time "0"
        tr_0 = pts.sel(time=timerange[0],method='nearest').dropna(dim="id").id.values
        tr_i = pts.sel(time=timerange[1],method='nearest').dropna(dim="id").id.values

        # get the difference between them (using sets)
        ids = list(set(tr_i)-set(tr_0))
        
        # apply the support selection function
        self._apply_selection(ids, 
                              selected_name=selected_name, 
                              replace_original=replace_original, 
                              selection_name=selection_name)
        return self
        
    
    def selectParticles_bycoords(self, xlim=None, zlim=None, tsel=None, 
                                 select_original=True,selected_name='',selection_name='',replace_original=False):
        '''
        coords : list = [[xmin, xmax],[zmin,zmax]]
        '''
        
        if xlim is None: xlim = self.xlimits
        if zlim is None: zlim = self.zlimits

        pts, selected_name = self._get_pts(select_original, selected_name)
            
        if tsel is None: tsel = 0
        
        pts = pts.sel(time=tsel, method='nearest')
        condX = ((pts["x"] >= xlim[0]) & (pts["x"] <= xlim[1])).compute()
        pts =  pts.where(condX, drop=True)
        
        condZ = ((pts["z"] >= zlim[0]) & (pts["z"] <= zlim[1])).compute()
        pts =  pts.where(condZ, drop=True)
        ids = pts.id.values
              
        # apply the support selection function
        self._apply_selection(ids, 
                              selected_name=selected_name, 
                              replace_original=replace_original, 
                              selection_name=selection_name)
        
        return self
    
    
    def selectParticles_bylayers(self, layers, tsel=None, 
                                select_original=True,selected_name='',selection_name='',replace_original=False):
        '''
        select particles by layer
        '''
        
        pts, selected_name = self._get_pts(select_original, selected_name)
        if tsel is None: tsel = 0 # Future: to create a function to evaluate an automatic tsel
        
        pts = pts.sel(time=tsel, method='nearest')
        cond = pts.layer.isin(layers).compute()
        ids =  pts.id.where(cond, drop=True).values
        
        self._apply_selection(ids, 
                              selected_name=selected_name, 
                              replace_original=replace_original, 
                              selection_name=selection_name)
        return self
    
    
    def classify_ParticlesRange(self, domain_intervals, tsel=None,
                               select_original=True,selected_name='',selection_name='',replace_original=False):
        
        #Classify all particles based on X ranges, given a time step tsel
        #Categories are based on the domain intervals keys
        
        pts, selected_name = self._get_pts(select_original, selected_name)
        
        if replace_original==True:
            selection_name = 'original'
        elif selection_name == '':
            selection_name == 'selected'
        
        if tsel is None: tsel = 0
        
        domain_intervals = domain_intervals.copy()
        
        try: field_name = domain_intervals['field_name']
        except: field_name = 'domain'
        del domain_intervals['field_name']
        
        typename = type(list(domain_intervals.keys())[1])
        if typename is str: typename='U256'
        
        if field_name in pts:
            print(f"{field_name} is a dataset variable, choose another name")
            return False
        
        pts[field_name] = (['id'], np.full(pts.sizes['id'], '', dtype=typename))
        
        snapshot = pts.sel(time=tsel, method='nearest')
        
        for dom, intervals in list(domain_intervals.items()):
            
            if not isinstance(intervals[0], (list, tuple, np.ndarray)): intervals = [intervals]
            mask_combined = np.zeros(snapshot.sizes['id'], dtype=bool)
            
            for start, end in intervals:
                mask_current = ((snapshot.x >= start) & (snapshot.x <= end)).compute()
                mask_combined |= mask_current.values
            
            ids_in_domain = snapshot.id.values[mask_combined]
            pts[field_name] = xr.where(pts.id.isin(ids_in_domain), dom,  pts[field_name])
        
        self.DTree.particles[selection_name] = pts
        
        pts[field_name].attrs['reference timestep'] = f'{tsel}myr'
        pts[field_name].attrs['classes range'] = str(domain_intervals)
        gc.collect()
        
        return self

    def fieldToParticle(self, variable, 
                       select_original=True,selected_name='',selection_name='',replace_original=False):
        
        pts0, selected_name = self._get_pts(select_original, selected_name)

        if replace_original==True:
            selection_name = 'original'
        elif selection_name == '':
            selection_name == 'selected'
        
        field = self.vars_DS[variable]
        #if not isinstance(variables, (list, tuple, np.ndarray)): variables = [variable]
        
        pts = pts0.copy()
        
        components = list(field.data_vars)
        
        field_interpolated = field.interp(
            x=pts.x, 
            z=pts.z, 
            method='linear'
        )
        
        for comp in components:
            pts[comp] = field_interpolated[comp].drop_vars(['x', 'z']).transpose('id', 'time') #did it worked?
        
        
        self.DTree.particles[selection_name] = pts

        gc.collect()
        return self

