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
from ._aux_functions import read_params

#Mandyoc Scenario class
class MandyocScen:
    """
    Core data manager for Mandyoc simulations (results).

    This class parses output directories, reads model parameters, and constructs 
    a unified DataTree. The classes alsos contains the scenarios params.txt 
    [future implementations will read temporal conditions]. 
    
    It aligns meshes (common variables), upscaled lithology, surface topography,
    and Lagrangian particles into a single navigable object. 

    Functions to manipulate data of `MandyocScen` are built-in and in acessors.

    Parameters
    ----------
    path : str
        The directory containing the Mandyoc output files in .nc and `param.txt`.
    
    variables : str or list of str, optional
        A list of standard variables to load (e.g., 'density', 'temperature'). 
        **Cannot be empty**. 
        
        Default is `['density']`.
    name : str, optional
        A custom name for the scenario. If None, the name of the base directory 
        is used. 
        
        Default is None.
    load_lithology : bool, optional
        If True, loads the upscaled lithology mesh into the `'/mesh/upscaled'` node. 
        
        Default is False.
    load_surface : bool, optional
        If True, loads 1D surface topography data into the `'/surface/topography'` 
        node. 
        
        Default is False.
    load_particles : bool, optional
        If True, loads Lagrangian particle tracking data into the `'/particles'` 
        node. 
        
        Default is False.
    particles_file : str, optional
        The filename of the NetCDF file containing particle trajectories. Only
        relevant if `load_particles` is True.
        
        Default is `'particles_trajectories.nc'`.
    xlimits : list of float, optional
        Spatial bounds [xmin, xmax] for data reading. If None, uses the 
        domain boundaries. Useful for huge scenarios.
        
        Default is None.
    zlimits : list of float, optional
        Spatial bounds [zmin, zmax] for data reading. If None, uses the 
        domain boundaries. Useful for huge scenarios.
        
        Default is None.
    tlimits : list of float, optional
        Temporal bounds [tmin, tmax] for data reading. If None, uses the 
        available time steps. Useful for long scenarios.
        
        Default is None.
    thick_air : float, optional
        The thickness of the sticky air layer (in meters). Default is 40e3.

    chunks_vars : dict, optional
        A dictionary defining the Dask chunking strategy for the spatial variables. 
        
        Default is `{"x": 'auto', "z": 'auto', 'time': "auto"}`.
    filter_air : bool, optional
        If True, filters out particles residing within the sticky air layer upon 
        loading. Only relevant if `load_particles` is True. 
        
        Default is True.
    air_layer : int, optional
        The numerical value or threshold identifying the air layer phase to be 
        filtered. Required if `filter_air` is True.
        
        Default is None.
    verbose : bool, optional
        If True, prints status messages to the console during loading and processing
        tasks. 
        
        Default is False.

    Attributes
    ----------
    path : pathlib.Path
        The absolute path to the scenario directory.
    name : str
        The designated name of the scenario.
    params : dict
        A dictionary of the parameters parsed from `param.txt`.
    DTree : xarray.DataTree
        The hierarchical data structure containing all loaded model outputs.
    xlimits : list of float
        The horizontal spatial limits applied to the data.
    zlimits : list of float
        The vertical spatial limits applied to the data.
    tlimits : list of float
        The temporal limits applied to the data.
    thick_air : float
        The sticky air thickness.
    z_corrected : bool
        Flag indicating whether the Z-coordinates have been corrected for topography.
    particles_loaded : bool
        Flag indicating whether the Lagrangian particles have been successfully loaded.
    
    Notes
    -----
    Check user guides for examples and tutorials. A comprehensive explanaition will be 
    available at Data Abstractions section.
    """

    def __init__(self, path:str, variables:list=['density'], name:str=None,
                 load_lithology:bool=False, load_surface:bool=False, load_particles:bool=False,
                 particles_file:str='particles_trajectories.nc',
                 xlimits:list=None, zlimits:list=None, tlimits:list=None, # ylimits should be implemented for the 3D version
                 thick_air:float=40e3,
                 chunks_vars:dict={"x": 'auto', "z": 'auto', 'time': "auto"},
                 filter_air:bool=True, air_layer:int=None, #only relevant if load particles
                 verbose:bool=False):

        # Setting directories and scen name
        self.path = Path(path)
        self.params = read_params(os.path.join(path,'param.txt'))
        self.verbose = verbose
        if isinstance(name, str): self.name = name
        else: self.name = self.path.name

        if self.verbose:
            print(f'Scenario at: {self.path}')
            print(f'Scenario name: {self.name}')
            print(f'Params - Ok')
        
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
        self.DTree = DataTree.from_dict( {"mesh":None,"surface":None,"particles":None})

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
            self.DTree['/mesh/original'] = xr.merge(standard_datasets, join='outer')

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
            self._load_particles(particles_file, chunks={'id': 'auto'}, filter_air=filter_air, air_layer=air_layer)
            
        if self.verbose: print(f"Particles [{particles_file}] loaded")
        
        return None
    
    def get_scenarioData(self, var:str='density'):
        """
        Extracts spatial and temporal metadata from a reference NetCDF file. 
        (Nx, Nz) are the grid dimension (elements) and 
        (XMAX, XMIN, ZMAX, ZMIN, TMAX, TMIN) are the dimensions boundaries.

        The extracted values are stored directly as class attributes.

        Parameters
        ----------
        var : str, optional
            The name of the variable to read (without the '.nc' extension). 
            
            Default is 'density'.

        Returns
        -------
        bool
            Returns True if the metadata was successfully extracted.

        Raises
        ------
        FileNotFoundError
            If the corresponding NetCDF file is not found in the scenario path.

        Notes
        -----
        This method is called automatically during class initialization to map 
        the scenario boundaries before building the full DataTree.
        """

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
        """
        Adjusts the Z-coordinates across the entire DataTree by a 
        uniform vertical shift to the Eulerian meshes, surface topography
        arrays, and particles based on the thickness of the sticky 
        air layer (or apply a custom baseline). 

        Useful to analyse scenario considering the baseline as 0 m.

        Parameters
        ----------
        factor : float, optional
            The value to subtract from the Z-coordinates (in meters). 
            If None, it is calculated automatically as
            `(self.ZMAX - self.thick_air)`. 
            
            Default is None.

        Returns
        -------
        bool
            Returns True if the correction was successfully applied. Returns 
            False if the coordinates were already corrected.

        Notes
        -----
        **Tree Mutation:** This method modifies the underlying datasets in-place. 
        It updates the `z` coordinates for nodes under `'/mesh'` and `'/particles'`, 
        and the `'surface'` data array for nodes under `'/surface'`. 
        
        **Attribute Mutation:** Updates `self.zlimits`, `self.ZMAX`, and 
        `self.ZMIN`, and sets `self.z_corrected` to True.
        """
        
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
        Internal method to load Eulerian grids and slice them to `xlimits`,`zlimits`.

        Parameters
        ----------
        variable : str
            The name of the variable to load.
        chunks: dict, optional
            Dictionary with the chunksizes (keys) for each dimensions (keys).

        Returns
        -------
        xr.Dataset
            The loaded data for the variable
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
        Loads Lagrangian particle trajectories from a NetCDF file into the DataTree.

        This internal method reads the particle file, applies lazy Dask 
        chunking, and truncates the dataset to match the active temporal limits 
        (`self.tlimits`). 
        It can also optionally filter out particles belonging to the sticky air
        phase before storing them in the tree.

        Parameters
        ----------
        name : str
            The filename of the particles NetCDF file to load.
        chunks : dict, optional
            A dictionary defining the Dask chunking strategy for the particle 
            dataset. 
            
            Default is `{'id': 'auto'}`.
        filter_air : bool, optional
            If True, evaluates the 'layer' variable to mask out particles of sticky air.
            
            Default is True.
        air_layer : int, optional
            The ID of the sticky air layer. If None and `filter_air` is True, the 
            algorithm assumes the maximum value in the 'layer' array corresponds
            to the air phase. 
            
            Default is None.

        Returns
        -------
        bool
            Returns True if the particles were successfully loaded, filtered, 
            and stored. Returns False if the specified file does not exist.

        Notes
        -----
        **Tree Mutation:** This method creates or overwrites the `'/particles/original'` 
        node in the `self.DTree` with the processed `xarray.Dataset`.
        
        **Attribute Mutation:** Sets the class attribute `self.particles_loaded` to True.
        """
        file_path = self.path / name
        
        if not file_path.exists():
            print(f"Warning: Particles file {file_path} not found. Skipping.")
            return False
            
        particles = xr.open_dataset(file_path, chunks=chunks)
        particles = particles.sel(time=slice(self.tlimits[0], self.tlimits[-1]))
        
        if filter_air and 'layer' in particles.data_vars:
            if isinstance(air_layer, int): air = air_layer
            else: air = int(particles.layer.max()) #future: find air layer by density?
                
            cond = particles.layer != air
            particles = particles.where(cond)
            if self.verbose: print(f'Air particles filtered [{air}]')
            
        self.DTree['/particles/original'] = particles
        self.particles_loaded = True
        return True
    
    # Future improvements: remove replace original, redudant if you can select "original"
    # Replace "selected" and "selection" for "source" and "target/dest"
    def _apply_selection(self, valid_ids, selected_name='', selection_name=''):
        """
        Internal function to select a dataset by the IDs and store it in the DataTree.
        """

        paths, names, pts = self._eval_particles_selection(selected_name=selected_name,selection_name=selection_name,get_selection=True)
        source_path, target_path = paths
        selected_name, selection_name = names

        if self.verbose:
            print(f'Selected particles: {selected_name}')
            print(f'Selection at: {selection_name}')

        # Lazy selection
        pts = self.DTree[source_path].ds
        self.DTree[target_path] = pts.sel(id=valid_ids)
    
        gc.collect()
        return None

    def _eval_particles_selection(self,selected_name, selection_name, get_selection=False):
        '''
        [Ongoing]: Evaluate the selection names and paths  -for particles
        '''

        if len(selection_name) == 0:
            selection_name = 'selected'
        
        elif len(selected_name) == 0:
            selected_name = 'original'
        
        source_path = f'/particles/{selected_name}'
        target_path = f'/particles/{selection_name}'

        if get_selection:
            return (source_path,target_path), (selected_name,selection_name), self.DTree[source_path].ds
        else:
            return (source_path,target_path), (selected_name,selection_name)

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
        paths, names,pts = self._eval_particles_selection(selected_name=selected_name,selection_name=selection_name,get_selection=True)
        source_path, target_path = paths
        selected_name, selection_name = names
        
        # select particles based on two snapshots
        snap_0 = pts.sel(time=timerange[0], method='nearest').compute()
        snap_i = pts.sel(time=timerange[1], method='nearest').compute()

        # selecting particles based on the X coordinate and discarding null values
        tr_0 = snap_0.id.where(~snap_0.x.isnull(), drop=True).values
        tr_i = snap_i.id.where(~snap_i.x.isnull(), drop=True).values

        # get the difference between them (using sets)
        ids = list(set(tr_i)-set(tr_0))

        if self.verbose:
            print(f'{len(ids)} particles were selected between {timerange[0]}-{timerange[1]}')
            
        # apply the support selection function
        self._apply_selection(ids, 
                              selected_name=selected_name,
                              selection_name=selection_name)
        return self
        
    
    def selectParticles_bycoords(self, xlim=None, zlim=None, tsel=None, 
                                 select_original=True,selected_name='',selection_name='',replace_original=False):
        '''
        coords : list = [[xmin, xmax],[zmin,zmax]]
        '''
        
        if xlim is None: xlim = self.xlimits
        if zlim is None: zlim = self.zlimits

        paths, names, pts = self._eval_particles_selection(selected_name=selected_name,selection_name=selection_name,get_selection=True)
        source_path, target_path = paths
        selected_name, selection_name = names
            
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
                              selection_name=selection_name)
        
        return self
    
    
    def selectParticles_bylayers(self, layers, tsel=None, selected_name='',selection_name=''):
        '''
        select particles by layer
        '''
        
        paths, names,pts0 = self._eval_particles_selection(selected_name=selected_name,selection_name=selection_name,get_selection=True)
        selected_name, selection_name = names
        pts = pts0.copy()

        if tsel is None: tsel = 0 # Future: to create a function to evaluate an automatic tsel
        
        pts = pts.sel(time=tsel, method='nearest')
        cond = pts.layer.isin(layers).compute()
        ids =  pts.id.where(cond, drop=True).values
        
        self._apply_selection(ids, 
                              selected_name=selected_name, 
                              selection_name=selection_name)
        return self
    
    
    def classify_ParticlesRange(self, domain_intervals, tsel=None,
                               selected_name='',selection_name=''):
        
        """
        [ongoing]
        """
        #Classify all particles based on X ranges, given a time step tsel
        #Categories are based on the domain intervals keys
        
        paths, names,pts0 = self._eval_particles_selection(selected_name=selected_name,selection_name=selection_name,get_selection=True)
        source_path, target_path = paths
        selected_name, selection_name = names
        pts = pts0.copy()
        
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
        
        target_path = f'/particles/{selection_name}'
        self.DTree[target_path] = pts
        
        pts[field_name].attrs['reference timestep'] = f'{tsel}myr'
        pts[field_name].attrs['classes range'] = str(domain_intervals)
        gc.collect()
        
        return self

    def fieldToParticle(self, variable, selected_name='', selection_name='', method='linear'):
        """
        [ongoing]
        """
        paths, names,pts0 = self._eval_particles_selection(selected_name=selected_name,selection_name=selection_name,get_selection=True)
        source_path, target_path = paths
        selected_name, selection_name = names
        pts = pts0.copy()

        field = self.DTree.mesh.original[variable]
        
        #if not isinstance(variables, (list, tuple, np.ndarray)): variables = [variable]
        
        #components = list(field.data_vars)
        
        field_interpolated = field.interp(
            x=pts.x, 
            z=pts.z,
            time=pts.time,
            method=method
        )
    
        pts[variable] = field_interpolated.drop_vars(['x', 'z']).transpose('id', 'time')

        self.DTree[source_path][variable] = pts

        gc.collect()
        return self
