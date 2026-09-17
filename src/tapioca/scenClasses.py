import os, gc, json, glob
from pathlib import Path
from shapely.geometry import Polygon,LineString
from shapely import to_ragged_array
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

import numpy as np
import xarray as xr
from xarray import DataTree

from ._variables import VARS_TYPES, VARIABLES_LIST, INTERFACES_PARAMETERS,DEFAULT_MATERIAL,MATERIAL_PARAMETERS,PARAMETERS_UNITS,SEC_PER_YEAR
from ._aux_functions import read_params, ensure_directory_exists, _numba_diffusion_loop

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
            The value to sum into the Z-coordinates (in meters).
            If None, it is calculated automatically as
            `(self.ZMAX + self.thick_air)`. 
            
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
            factor = self.thick_air 
        
        # Traverse specific tree nodes and apply Z-correction if 'z' exists
        structs_to_check = ['/mesh', '/surface','/particles']
        
        mesh_nodes = list(self.DTree['/mesh'].children) # Mesh data (Eulerian grids)
        surface_nodes = list(self.DTree['/surface'].children) # Surface/string data (only X dimension)
        particles_nodes = list(self.DTree['/particles'].children) # Lagrangian particles 
        
        for node in mesh_nodes:
            node_ds = self.DTree['/mesh'][node].ds
            self.DTree['/mesh'][node] = node_ds.assign_coords(z=(node_ds['z'] + factor))
        
        # Add a If condition for when the used does not load all data types
        for node in surface_nodes:    #O script em julia está exportando com a superfície corrigida, mudar para exportar com o dado ORIGINAL
            node_ds = self.DTree['/surface'][node].ds
            self.DTree['/surface'][node] = node_ds.assign(surface=(node_ds['surface'] + factor))

        for node in particles_nodes:
            node_ds = self.DTree['/particles'][node].ds
            self.DTree['/particles'][node] = node_ds.assign(z=(node_ds['z'] + factor))

        # The dataset attributes must be updated too
        self.zlimits = [self.zlimits[0] + factor, self.zlimits[1] + factor]
        self.ZMAX += factor
        self.ZMIN += factor
        
        if self.verbose:
            print(f"Z coordinate corrected by + {factor} m")
            print(f"New z limits: {self.zlimits}")
            
        self.z_corrected = True
        return True

    def filter_air_fields(self, dens_threshold=10.0, air_value=-1):
        """
        Filters variable fields of the sticky air layer by density.
        All densities below a threshold will be considered air. 

        Parameters
        ----------
        dens_threshold : float
            The maximum density to consider the air layer, in kg/m^3.

            Default is 10.
        air_value : int or float
            The new value into air layer.

            Default is -1.
        """
        return

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


class MandyocBuilder:

    def __init__(self, path:str, Nx:int, Nz:int, Lx:float, Lz:float, 
            Ny:int=0, Ly:float=0.0, thick_air:float=40e3, 
            verbose:bool=True):
        """
        A builder class made to generate, manage, and export initial setups 
        for Mandyoc geodynamic models. 
        
        It utilizes xarray's DataTree structure to strictly separate 
        1D geometric boundaries (interfaces) and 2D/3D material properties (fields). 
        By default, it assumes a Cartesian coordinate system where the Z-axis is negative 
        downwards (from the bottom of the model at -Lz up to 0). 

        Parameters
        ----------
        path : str or Path
            The absolute or relative directory path where the generated scenario files will be exported.
        
        Nx : int
            The number of numerical nodes along the X-axis (horizontal).
            
        Nz : int
            The number of numerical nodes along the Z-axis (vertical/depth).
            
        Lx : float
            The total length of the model domain along the X-axis (in meters).
            
        Lz : float
            The total depth of the model domain (in meters). 
            Note: The class automatically converts this to a negative value internally to represent depth below the surface.
            
        Ny : int, optional
            The number of numerical nodes along the Y-axis. Required only for 3D scenarios. 
            Default is 0.
            
        Ly : float, optional
            The total physical length of the model domain along the Y-axis in meters. 
            Required only for 3D scenarios. 
            Default is 0.0.
            
        thick_air : float, optional
            The thickness of the sticky air layer at the top of the model domain (in meters). 
            Default is 40000.0 (40 km).
            
        verbose : bool, optional
            If True, the builder will print status updates and configuration details to the console during functions execution. 
            Default is True.

        Attributes
        ----------
        path : str
            The output directory path.
            
        verbose : bool
            The current verbosity state of the builder.
            
        thick_air : float
            The defined thickness of the sticky air layer.
            
        Nx, Nz, Ny : int
            The node counts for the X, Z, and Y axes, respectively. 
            `Ny` remains None if the model is 2D.
            
        Lx, Lz, Ly : float
            The physical dimensions of the grid. `Lz` is strictly negative. 
            `Ly` remains None if the model is 2D.
            
        dimensions : int
            The spatial dimensionality of the model (either 2 or 3), determined automatically 
            based on the provided `Ny` and `Ly` parameters.
            
        x : numpy.ndarray
            A 1D array of length `Nx` containing the spatial coordinates along the X-axis 
            (from 0 to `Lx`).
            
        z : numpy.ndarray
            A 1D array of length `Nz` containing the spatial coordinates along the Z-axis 
            (from `Lz` to 0).
            
        y : numpy.ndarray or None
            A 1D array of length `Ny` containing the spatial coordinates along the Y-axis 
            (from 0 to `Ly`). None if the model is 2D.
            
        DTree : xarray.DataTree
            The core hierarchical data structure storing the model's geometry and physics. 
            It contains two main groups:
            * `/interfaces`: An `xarray.Dataset` storing the 1D (or 2D if the model is 3D) 
            depth arrays defining the boundaries between different geological layers. 
            Contains material parameters in the datarrays attributes.
            * `/fields`: An `xarray.Dataset` storing the fully 2D (or 3D) grids containing 
            the final material properties (e.g., density, radiogenic heat) in the discrete
            numerical nodes. It can also contain the input variables (vx, vz, temperature).
            """
        
        self.path = path
        self.verbose = verbose
        self.thick_air = thick_air
        ensure_directory_exists(self.path)

        self.Nx = Nx
        self.Nz = Nz
        self.Lx = Lx
        self.Lz = -Lz
        self.Ly = self.Ny = self.y = None
        self.dimensions = 2

        self.x = np.linspace(0, self.Lx, self.Nx)
        self.z = np.linspace(self.Lz, 0, self.Nz)

        # Creating the DataTree structure to store the model layers
        self.DTree = DataTree.from_dict( {"fields":
                            xr.Dataset(coords={'x': self.x, 'z': self.z}), 
                                    "interfaces":xr.Dataset(coords={'x': self.x})})

        self.DTree['/fields']["x"].attrs["units"] = "m"
        self.DTree['/fields']["z"].attrs["units"] = "m"
        self.DTree['/fields']["x"].attrs["axis"] = "X"
        self.DTree['/fields']["z"].attrs["axis"] = "Z"
        self.DTree['/fields'].attrs["description"] = "This group contains the fields of the model, each variable is a DataArray considering the material properties from the interfaces."

        self.DTree['/interfaces']["x"].attrs["units"] = "m"
        self.DTree['/interfaces'].attrs["description"] = "This group contains the interfaces of the model, each interface is a DataArray with the depth and the material properties as attributes."
        
        # setting the model dimensions
        if Ny>0 and abs(Ly)>0:
            self.Ly = Ly
            self.Ny = Ny
            self.dimensions = 3
            self.y = np.linspace(0, self.Ly, self.Ny)
            self.DTree = DataTree.from_dict( {"fields":
                                xr.Dataset(coords={'x': self.x, 'y': self.y, 'z': self.z}), 
                                        "interfaces":xr.Dataset(coords={'x': self.x, 'y': self.y})})

            self.DTree['/fields'].expand_dims(dim='y', axis=1).assign_coords(y=self.y)
            self.DTree['/interfaces'].expand_dims(dim='y', axis=1).assign_coords(y=self.y)
            
            self.DTree['/fields']["y"].attrs["units"] = "m"
            self.DTree['/fields']["y"].attrs["axis"] = "Y"
            self.DTree['/interfaces']["y"].attrs["units"] = "m"
            
    def _print_verbose(self, message:str, ending ="\n"):
        '''
        Internal method to print a message if verbose is True

        Parameters
        ----------
        message : str
            The printed message

        ending : str
            The end of the string text.
            Default is `\n`.
        '''
        if self.verbose:
            print(message,end=ending)
        return None

    def create_interface(self, id_layer:int, position:float|np.ndarray, interface_name:str=None, **material):
        '''
        Create a structural interface and assign its geometric position and material properties.

        This function creates an interface DataArray and adds it to the model's DataTree. It supports both 2D and 3D dimensions. 
        Note that an ID of -1 is treated as a special uppermost bounding layer (e.g., air); its coordinates are forced to 0 regardless of the position input.
        
        Parameters
        ----------
        id_layer: int
            The integer identifier for the layer. Use -1 to specify the uppermost layer of the model.
        
        position: float or np.ndarray
            The depth position of the interface. This can be a constant float representing a flat depth, or an array defining a complex surface.
        
        interface_name: str, optional
            The name assigned to the interface in the dataset. If not provided, it automatically defaults to `interface_{id_layer}`.

            Default is None.

        **material
            Rheological and thermal parameters for the layer. These inputs overwrite the properties in DEFAULT_MATERIAL. 
            A density value ('rho') must be explicitly provided or the function will raise a ValueError.
        
        '''

        if interface_name==None:
            interface_name = f'interface_{id_layer}'
                
        # Defining material properties
        interface_parameters = DEFAULT_MATERIAL.copy() # default material

        self._print_verbose(f'ID {id_layer} -> {interface_name}')
        for k in list(material.keys()):
            interface_parameters[k] = material[k]

            self._print_verbose(f'{k}:{interface_parameters[k]}', ending="; ")
        self._print_verbose('')

        if interface_parameters['rho'] == None:
            raise ValueError("ERROR: You must to set a density (rho) to the material interface!")


        if id_layer!=-1:
            interface_coords = self._evaluate_interface_position(position)
            #create a default interface at the top of the model (ZMAX)

        else:
            print(f'''Layer with ID {id_layer} is the top layer, it contains the initial material properties of the model that is above the last interface.\nYou can set material properties for this layer, but you cannot create an interface.''')
            
            interface_coords = np.ones(self.Nx) if self.dimensions==2 else np.ones((self.Nx,self.Ny))
            interface_coords *= 0

        # Check interface_coords shape and values to ensure they are within the model bounds
        # self._check_coords(interface_coords)

        if self.dimensions==2:
            self.DTree['/interfaces'][interface_name] = (xr.DataArray(interface_coords, coords=[self.x], 
                                                                              dims=['x'], name=interface_name))
            
        elif self.dimensions==3:
            self.DTree['/interfaces'][interface_name] = (xr.DataArray(interface_coords, coords=[self.x,self.y], 
                                                                                 dims=['x','y'], name=interface_name))

        self.DTree['/interfaces'][interface_name].attrs['id'] = id_layer
        self.DTree['/interfaces'][interface_name].attrs.update(interface_parameters)
        return self


    def add_interface_from_shapely(self, layer:str, polygon:Polygon, poly_name:str=None, intersection_line:str='mid', **material): # 2D
        '''
        Add interfaces from a defined convex polygon. The user must specify the layer in which the polygon will be contained.

        This function creates the interfaces DataArray and adds it to the model's DataTree. 
        Note that this function creates two interfaces: top and bottom of the polygon. The IDs are automatically corrected using the specified layer.
        The top interface receives the material within the polygon, while the bottom interface receives the parameters from the outer layer.
        
        Parameters
        ----------
        layer: str
            The name of the layer that will receive/contain the given polygon.
        
        polygon: shapely.Polygon
            The Polygon object created with shapely. The polygon has to be convex and with no negative features in x direction. 
        
        poly_name: str, optional
            The name for the polygon interfaces.

            Default is None.
        
        intersection_line: str, optional
            The relative position of the interface line in the polygon. Possible values are 'mid', 'top', 'bot'.

            Default is 'mid'.

        **material
            Rheological and thermal parameters within the polygon.
        
        '''
        polycoords = to_ragged_array([polygon])[1]
        polyxcoords = np.unique(polycoords[:,0])

        xmin,zmin, xmax,zmax = polygon.bounds
        mid_z = (zmin+zmax)/2

        # checking where the null line position
        if intersection_line == 'top':
            z_ref = zmax
        elif intersection_line == 'bot':
            z_ref = zmin
        elif intersection_line == 'mid':
            z_ref = mid_z
        else:
            raise ValueError("The position for the intersection line must be: `top`, `bot` or `mid`.")

        self._print_verbose(f'Polygon bounds: ({xmin,zmin}) ({xmax,zmax})')
        self._print_verbose(f'Intersection line ({intersection_line} line) is at {z_ref} m depth')

        top_intersec_pts = []
        bot_intersec_pts = []

        for x in polyxcoords:
            vline = LineString([(x,0),(x,self.Lz)])
            pts_x, pts_y = vline.intersection(polygon).xy
            pts_x, pts_y = list(pts_x),list(pts_y)
            if len(pts_x) > 1:
                top_intersec_pts.append([pts_x[0],pts_y[0]]) # top point
                bot_intersec_pts.append([pts_x[1],pts_y[1]]) # bot point

            elif len(pts_x) == 1: # append the same point 
                top_intersec_pts.append([pts_x[0],pts_y[0]])
                bot_intersec_pts.append([pts_x[0],pts_y[0]])

        fac = 1 # m
        special_points = [ 
        [xmin-fac, z_ref], # left polygon boundary
        [xmax+fac, z_ref], # right polygon boundary
        [0, z_ref], # left boundary
        [self.Lx, z_ref], # right boundary
        ]

        for sp in special_points:
            top_intersec_pts.append(sp)
            bot_intersec_pts.append(sp)

        top_interface = np.array(top_intersec_pts)
        bot_interface = np.array(bot_intersec_pts)
        
        # handle with the ids
        interfaces_names = list(self.DTree.interfaces.variables)[:-1]
        id_above = self.DTree['/interfaces'][layer].attrs['id']

        for n in interfaces_names:
            if self.DTree['/interfaces'][n].attrs['id'] == id_above:
                above_material = self.DTree['/interfaces'][n].attrs.copy()
                self._print_verbose(f"Above layer is: ({self.DTree['/interfaces'][n].attrs['id']}) {n}")

            if (self.DTree['/interfaces'][n].attrs['id'] >= id_above) and (id_above > 0):
                self.DTree['/interfaces'][n].attrs['id'] += 2

        # in case the polygon is within the air (-1)
        if id_above < 0:
            id_above = len(interfaces_names)+1
        
        self.create_interface(id_above,bot_interface,interface_name=f'bot_{poly_name}',**above_material)
        self.create_interface(id_above+1,top_interface,interface_name=f'top_{poly_name}',**material)

        return self

    def _evaluate_interface_position(self, coord):
        """
        Internal method to verify what type of coordinate or position the interfaces function
        is recieving
        """

        if isinstance(coord,float) or isinstance(coord,int):
            self._print_verbose('The interface is a flat line/surface')
        
            coordArray = np.ones(self.Nx) if self.dimensions==2 else np.ones((self.Nx,self.Ny))
            return coordArray*coord

        # Verify if the coord is a bunch of points to interpolate or is the whole interface (line/surface)
        elif isinstance(coord,np.ndarray) and self.dimensions==2:
            self._print_verbose('The interface is a line')

            if len(coord) == len(self.x):
                print('The whole interfaces was provided, no need to interpolate')
                #print(f"WARNING: The number of points provided ({len(coord)}) is different from the number of x-coordinates ({self.Nx}).)
                return coord
            
            elif len(coord[0])!=2:
                raise ValueError("ERROR: The coordinates for the interface must be a 2D array with two columns (x,z).")
            args_sorted = np.argsort(coord[:,0])

            return np.interp(self.x, coord[:,0][args_sorted], coord[:,1][args_sorted])

        else:
            raise ValueError("ERROR: The coordinates for the interface must be: " \
            "float;" \
            "int for 2D and 3D models; " \
            "2D array with points (x,z columns) for 2D models.")

        # Add interfaces from Shapely geometries
        
        return None


    def evaluate_interfaces(self):
        #analyse if interfaces are crossing each other and to fix it

        return self

    def create_materials_fields(self):
        """
        Turns interfaces into 2D (or 3D) material fields.

        This method iterates through all defined material parameters (e.g., density, radiogenic heat) 
        and constructs corresponding spatial fields in the `DataTree`.

        The interfaces are automatically sorted by their `id` attribute before field generation. 
        The method handles two special boundary IDs to bound the top and bottom of the model:
        * ID `-1`: Represents the uppermost layer (e.g., sticky air). It fills all grid nodes 
          above the highest interface.
        * ID `0`: Represents the lowermost basal layer. It fills all grid nodes below the deepest interface.
        """
        #create the xdataarray fields with interfaces params 

        self._print_verbose("Sorting interfaces by ID:")
        sorted_interfaces = sorted(self.DTree.interfaces.ds.data_vars, key=lambda v: self.DTree.interfaces.ds[v].attrs['id'])
        self.DTree.interfaces.ds = self.DTree.interfaces.ds[sorted_interfaces]

        interfaces_names = list(self.DTree.interfaces.variables)[:-1]
        self._print_verbose(f"{'-'.join(interfaces_names)}")

        for p in MATERIAL_PARAMETERS:
            field = xr.DataArray(
                                np.ones((self.Nx, self.Nz))*-999, 
                                dims=('x', 'z'),
                                coords={ 'x':self.x, 'z': self.z}
                                )

            for i in range(len(interfaces_names)):
                name = interfaces_names[i]
                curr_interface = self.DTree.interfaces.ds[name]
                id_layer = curr_interface.attrs['id']
                fill_value = curr_interface.attrs[p]

                if id_layer < 1:
                    continue
                
                below_interface = self.DTree.interfaces.ds[interfaces_names[i-1]]
                id_layer_below = below_interface.attrs['id']

                cond = (field.z > below_interface) & (field.z <= curr_interface)

                field = xr.where(cond, fill_value, field)

            # Id 0
            name = interfaces_names[1]
            curr_interface = self.DTree.interfaces.ds[name]
            id_layer = curr_interface.attrs['id']
            fill_value = curr_interface.attrs[p]

            cond = (field.z <= curr_interface)
            field = xr.where(cond, fill_value, field)

            # Id -1
            name = interfaces_names[0]
            curr_interface = self.DTree.interfaces.ds[name]
            id_layer = curr_interface.attrs['id']
            fill_value = curr_interface.attrs[p]

            below_interface = self.DTree.interfaces.ds[interfaces_names[-1]]

            cond = (field.z > below_interface)
            field = xr.where(cond, fill_value, field)

            field.attrs['parameter'] = INTERFACES_PARAMETERS[p]
            field.attrs['unit'] = PARAMETERS_UNITS[p]

            self.DTree.fields[p] = field
            self._print_verbose(f'Field created: {INTERFACES_PARAMETERS[p]} [{PARAMETERS_UNITS[p]}]')

        return self

    def create_velocity_field(self, velocbuilder:VelocityFieldBuilder=None, 
                                  vxconst:float=0,vzconst:float=0):
        '''
        Gets the `VelocityFieldBuilder` class to create the vx and vz fields within the 
        scenario DataTree. Constant values of `vx` and `vz` can be given. 
        This function automatically converts cm/y to m/s.

        Parameters
        ----------
        velocbuilder:VelocityFieldBuilder
            The previous setup of the velocity field in the boundaries. It must be 
            already treated to avoid divergences.
        
        vx:float, optional
            Value for the vx if velocbuilder is not gave.

            Default is 0.
        
        vz:float, optional
            Value for the vz if velocbuilder is not gave.

            Default is 0.
        '''
    
        vx = xr.DataArray(
                        np.ones((self.Nx, self.Nz))*vxconst, 
                        dims=('x', 'z'),
                        coords={ 'x':self.x, 'z': self.z}
                        )
        
        vz = xr.DataArray(
                        np.ones((self.Nx, self.Nz))*vzconst, 
                        dims=('x', 'z'),
                        coords={ 'x':self.x, 'z': self.z}
                        )

        fac = 1
        if velocbuilder.units == 'cm/y': 
            fac = 1 / (100 * SEC_PER_YEAR)
            self._print_verbose('Velocity field are in cm/y, applying corrections...')

        if velocbuilder is not None:
            vx = xr.where(vx.z==0, velocbuilder.velocs.top.vx, vx).transpose('x', 'z')
            vx = xr.where(vx.z==self.Lz, velocbuilder.velocs.bot.vx, vx).transpose('x', 'z')
            vx = xr.where(vx.x==0, velocbuilder.velocs.left.vx, vx).transpose('x', 'z')
            vx = xr.where(vx.x==self.Lx, velocbuilder.velocs.right.vx, vx).transpose('x', 'z')
            
            vz = xr.where(vz.x==0, velocbuilder.velocs.left.vz, vz).transpose('x', 'z')
            vz = xr.where(vz.x==self.Lx, velocbuilder.velocs.right.vz, vz).transpose('x', 'z')
            vz = xr.where(vz.z==0, velocbuilder.velocs.top.vz, vz).transpose('x', 'z')
            vz = xr.where(vz.z==self.Lz, velocbuilder.velocs.bot.vz, vz).transpose('x', 'z')

            self.DTree.fields['vx'] = vx * fac
            self.DTree.fields['vz'] = vz * fac


        self.DTree.fields['vx'].attrs['unit'] = 'm/s'
        self.DTree.fields['vz'].attrs['unit'] = 'm/s'

        self._print_verbose('Velocity field was created in the scenario builder.')

        return self

    def create_temperature_field(self, tempbuilder:TemperatureFieldBuilder):
        '''
        Gets the `TemperatureFieldBuilder` class to create the temperature field within 
        the scenario DataTree.

        Parameters
        ----------
        tempbuilder:TemperatureFieldBuilder
            The previous setup of the temperature field. It must be already treated.
        '''
    
        temp = xr.DataArray(
                        np.ones((self.Nx, self.Nz))*-1, 
                        dims=('x', 'z'),
                        coords={ 'x':self.x, 'z': self.z}
                        )


        temp[:] = tempbuilder.temperature
        
        self.DTree.fields['temperature'] = temp
        self.DTree.fields['temperature'].attrs['unit'] = 'deg C'
        self.DTree.fields['temperature'].attrs['t_potential'] = tempbuilder.t_pot

        self._print_verbose('Temperature field was created in the scenario builder.')

        return self

    def export_interfaces(self, export:str='dataset'):
        '''
        Export the created interfaces into the "interfaces.txt" file required in Mandyoc.

        Parameters
        ----------

        export : str, optional
            Mode of saving the created interfaces into a more readable file. Support 
            modes are 'CSV' and 'dataset' (netcdf file).

            Default is `'dataset'`. 
        '''

        texts = {p: [] for p in MATERIAL_PARAMETERS}
        dic_params = {p: [] for p in MATERIAL_PARAMETERS}
        interfaces_list = []

        interfaces_names = list(self.DTree.interfaces.ds.data_vars)
        interfaces_names = interfaces_names[1:] + [interfaces_names[0]]
        
        self._print_verbose("Exporting dataset to interfaces.txt")
        self._print_verbose(f"Interfaces: {' -> '.join(interfaces_names)}")

        for name in interfaces_names:    
            layer = self.DTree.interfaces.ds[name]
            interface = layer.values
            material = layer.attrs
            for p in MATERIAL_PARAMETERS:
                texts[p].append(str(material[p]))

                if export.lower()=='csv':
                    dic_params[p].append(material[p])

            interfaces_list.append(interface)

        interfaces_list = np.array(interfaces_list).T

        base_text=f'''C  {"  ".join(texts['C'])}
        rho  {"  ".join(texts['rho'])}
        H  {"  ".join(texts['H'])}
        A  {"  ".join(texts['A'])}
        n  {"  ".join(texts['n'])}
        Q  {"  ".join(texts['Q'])}
        V  {"  ".join(texts['V'])}
        k  {"  ".join(texts['k'])}
        weakening_seed  {"  ".join(texts['weakening_seed'])}
        cohesion_min  {"  ".join(texts['cohesion_min'])}
        cohesion_max  {"  ".join(texts['cohesion_max'])}
        friction_angle_min  {"  ".join(texts['friction_angle_min'])}
        friction_angle_max  {"  ".join(texts['friction_angle_max'])}
        '''

        max_seq = 0
        for line in base_text.split("\n"):
            for seq in line.strip().split('  ')[1:]:
                # print(seq)
                max_seq = len(seq.strip()) if len(seq.strip()) > max_seq else max_seq

        max_seq += 3
        print(max_seq)
        with open(f"{self.path}interfaces.txt", "w") as f:
            for line in base_text.split("\n")[:-1]:
                all_seqs = line.strip().split('  ')
                f.write(all_seqs[0].ljust(21)) # 21 spaces -> len of "friction_angle_***   "

                for seq in all_seqs[1:]: # add parameters text
                    f.write(''.join(seq).ljust(max_seq))

                f.write("\n")

            np.savetxt(f, interfaces_list, fmt="%.1f")

            f.close()

        if export.lower() == 'dataset':
            self._print_verbose('Exporting interfaces to NETCDF file')
            nc_path = os.path.join(self.path, "interfaces_parameters.nc")
            self.DTree.interfaces.ds.to_netcdf(nc_path)
    
        elif export.lower() == 'csv':
            from pandas import DataFrame
            self._print_verbose('Exporting interfaces to CSV file')
            dic_params['name'] = interfaces_names
            df = DataFrame(dic_params)
            csv_path = os.path.join(self.path, "interfaces_parameters.csv")
            df.set_index('name').T.to_csv(csv_path, sep=';')
            
        return True

    def export_field(self,field:str,header:str=''):
        '''
        Function to export fields in the mandyoc required format.
        For the velocity or temperature field, this function exports the field with the 
        appropriates name.

        Parameters
        ----------
        field:str
            Field to be exported. The special fields are (`velocity`,`temperature`)
        header:str, optional
            Comments in the header of the file. Default is ''.
        '''
        name=field
        if len(header)==0: header='v1\nv2\nv3\nv4'
        
        if field == 'velocity':    
            vx = self.DTree.fields.vx
            vz = self.DTree.fields.vz

            vvx = vx.values.reshape(self.Nx*self.Nz)
            vvz = vz.values.reshape(self.Nx*self.Nz)

            velocity_export = np.zeros((2, self.Nx * self.Nz))
            velocity_export[0,:] = vvx
            velocity_export[1,:] = vvz

            data_export = np.reshape(velocity_export.T, (np.size(velocity_export)))
            name = 'input_velocity_0'

        else:
            data = self.DTree.fields[field]
            data_export = np.reshape(data.values, (self.Nx * self.Nz))
            
            if field=='temperature':
                name='input_temperature_0'

        self._print_verbose(f'Exporting velocity field ({len(data_export)})')
        np.savetxt(f"{self.path}{name}.txt", data_export, header=header)

        return self
    
    def set_params(self, params_str:str='', **kwargs):
        """
        Set parameters for the model
        """

        # create a function to evaluate params!
        self.params = kwargs
        
        return self

    def plot_interfaces(self, fig:plt.Figure, axes:plt.Axes, mode:str='fill'):

        return fig, axes



class VelocityFieldBuilder:

    def __init__(self, scenarioBuilder:MandyocBuilder, units='cm/y'):
        '''
        Class to build a correctly and conservative velocity field for a mandyoc scenario.
        
        The basic workflow of this class is:
        (1) Initialize the VelocityFieldBuilder with your current scenario;
        (2) Set velocities using `set_region` or `linear_velocity` functions;
        (3) Check the field conservation using the function `integrate_normal_components`;
        (4) Apply the velocity correction using the function `conservate_veloc`. Maybe you 
        should to iterate this step until to reach a tolerance divergence.
        (5) Brief visualize the velocities with `view_veloc`.

        Parameters
        ----------

        scenarioBuilder: MandyocBuilder
            The scenario that will recieve the velocity field. The VelocityFieldBuilder will 
            get all parameters and attributes from the scenario class.
        
        units: str, optional
            The units of the velocity field. It must be 'cm/y' or 'm/s'.

            Default is 'cm/y'.

        Attributes
        ----------
        scenario: MandyocBuilder
        Nx: int
        Nz: int
        Lx: float
        Lz: float
        x: numpy.ndarray
        z: numpy.ndarray
        boundaries: list
        velocs: xarray.DataTree
        units: str

        Ny: int
        Ly: float
        y: numpy.ndarray

        '''
        self.scenario = scenarioBuilder

        self.Nx = scenarioBuilder.Nx
        self.Nz = scenarioBuilder.Nz
        self.Lx = scenarioBuilder.Lx
        self.Lz = scenarioBuilder.Lz

        self.x = scenarioBuilder.x
        self.z = scenarioBuilder.z

        self.dx = self.x[1] - self.x[0]
        self.dz = self.z[1] - self.z[0]

        self.boundaries = ['left','right','top','bot']
        #self.velocs = xr.Dataset(dims=('x','z'),coords={'x':self.x, 'z':self.z})

        baseZ_array = xr.DataArray((np.zeros(self.Nz)),dims=('z'),coords={'z':self.z})
        baseX_array = xr.DataArray((np.zeros(self.Nx)),dims=('x'),coords={'x':self.x})

        self.velocs = DataTree.from_dict( {'left':xr.Dataset({'vx':baseZ_array,'vz':baseZ_array},coords={'z':self.z}),
                             'right':xr.Dataset({'vx':baseZ_array,'vz':baseZ_array},coords={'z':self.z}),
                             'top':xr.Dataset({'vx':baseX_array,'vz':baseX_array},coords={'x':self.x}),
                             'bot':xr.Dataset({'vx':baseX_array,'vz':baseX_array},coords={'x':self.x})}
                             )
        

        if scenarioBuilder.dimensions == 3:
            self.Ny = scenarioBuilder.Ny
            self.Ly = scenarioBuilder.Ly 
            self.y = scenarioBuilder.y
            self.dy = self.y[1] - self.y[0]

            self.boundaries.append('front')
            self.boundaries.append('back')

            # The 3D still have to be implemented

        self.units = units
        scenarioBuilder._print_verbose(f'The velocity builder was created. Velocities have to be in {self.units}')

    def set_region(self, boundary:str, range:list|tuple|np.ndarray, vx=None, vz=None, mode:str='set'):
        '''
        This function set or add a constant velocity in a region of the selected boundary.

        Parameters
        ----------
        boundary: str
            The boundary that you recieve the imposed velocity ('top','bot','left','right).

        range: list|tuple|np.ndarray
            The min and max values on the given boundary, such as [XMIN,XMAX] or [ZMIN,ZMAX].

        vx: float, optional
            The velocity in the Vx component.

            Default is None.

        vz: float, optional
            The velocity in the Vz component.

            Default is None.
        
        mode: str, optional
            The mode that the velocity will be applied: 'add' or 'set'.
            
            Default is 'set'.
        '''
        mask = self._create_mask(values_range=range, boundary=boundary)

        self._apply_veloc(mask, boundary, vx, vz, mode)

        if vx is None: vx = 0
        if vx is None: vz = 0
        self.scenario._print_verbose(f'Constant velocity in the {boundary}: from [{range[0]:.1f},{range[1]:.1f}], vx={vx:.3e}, vz={vz:.3e} {self.units}')

        return self

    def linear_velocity(self, boundary:str, v_component:str, range:list|tuple|np.ndarray=None, 
                        clims:list|tuple|np.ndarray=None, vlims:list|tuple|np.ndarray=None, 
                        coef:float=None, b:float=None, mode:str='set'):
        '''
        This function apply a linear velocity function in a range of the selected boundary.
        
        V = a * COORD + b

        The user can give the limits and the velocities to calculate `a` and `b` or give 
        these parameters directly.  

        Parameters
        ----------
        boundary: str
            The boundary that you recieve the imposed velocity ('top','bot','left','right).

        v_component: str
            Which velocity component are being set ('vx' or 'vz').
            
        range: list|tuple|numpy.ndarray
            The min and max values on the given boundary to apply the velocities, such 
            as [XMIN,XMAX] or [ZMIN,ZMAX].

        clims: list|tuple|numpy.ndarray, optional
            The min and max values in the coordinate to calculate the line slope, such 
            as [XMIN,XMAX] or [ZMIN,ZMAX]. If None, `coef` and `b` must be given. 

            Default is None.

        vlims: list|tuple|numpy.ndarray, optional
            The velocity values in the coordinates given in `clims`. 
            If None, `coef` and `b` must be given. 

            Default is None.
        
        coef: float, optional
            The slope of the function. 

            Default is None.

        b: float, optional
            The constant of the function. 

            Default is None.

        mode: str, optional
            The mode that the velocity will be applied: 'add' or 'set'.
            
            Default is 'set'.
        '''

        vx = self.velocs[boundary].vx
        vz = self.velocs[boundary].vz 

        if coef is None: # calculate the slope coeficient
            coef = (vlims[1]-vlims[0])/(clims[1]-clims[0])
            
            if b is None: # the linear constant
                b = vlims[1] - coef * clims[1]

        else: # user gives the coeficient 
            if b is None: b = 0.0


        varying_coord = self._evaluate_bound(boundary)
        veloc = coef * varying_coord + b

        self.scenario._print_verbose(f'linear velocity in the {boundary} boundary: {v_component} = {coef:.3e} * COORD + {b:.3e}')

        if v_component == 'vz': 
            vz = veloc
        
        elif v_component == 'vx': 
            vx = veloc
        
        mask = self._create_mask(values_range=range, boundary=boundary)

        self._apply_veloc(mask, boundary, vx, vz, mode)
        
        return self

    def _create_mask(self, values_range, boundary):
        ''' 
        Internal method to get the masked interval
        '''

        if boundary in ['top','bot']:
            xmin, xmax = values_range
            mask = (self.velocs[boundary].x >= xmin) & (self.velocs[boundary].x <= xmax)

        elif boundary in ['left','right']:
            zmin,zmax = values_range
            mask = (self.velocs[boundary].z >= zmin) & (self.velocs[boundary].z <= zmax)

        return mask

    def integrate_normal_components(self):
        '''
        Gives outflux integration (of the normal components) of the velocity field by:

            Int V•n dS = [ Int Vx(x=Lx) dz - Int Vx(x=0) dz ] + [Int Vz(0) dx - Int Vz(Lz) dx ]

        returns: float
            Return the outflux.
        '''

        left_vx = self.velocs['left'].vx.sortby('z')
        right_vx = self.velocs['right'].vx.sortby('z')
        top_vz = self.velocs['top'].vz.sortby('x')
        bot_vz = self.velocs['bot'].vz.sortby('x')

        outflux_l = -left_vx.integrate("z")
        outflux_r = right_vx.integrate("z")

        outflux_t = top_vz.integrate("x")
        outflux_b = -bot_vz.integrate("x")

        sum_outflux = outflux_l + outflux_r + outflux_t + outflux_b
        self.scenario._print_verbose(f"Int V = {sum_outflux.values:.3e}")
        
        return sum_outflux

    def conservate_veloc(self,boundary:str='bot', correction_factor:float=1.0):
        '''
        Apply the divergence of the velocity field, or a fraction of it, field into 
        a whole boundary. 

        Parameters
        ----------
        boundary: str
            The boundary that you recieve the imposed velocity ('top','bot','left','right).
        
        correction_factor: float
            The fraction of correction to be applied. Should be: 0.0 < CF <= 1.0.

            Default is 1.0.
        '''

        v_exc = self.integrate_normal_components()
        target_flux = -v_exc * correction_factor

        if boundary in ['bot', 'top']:
            normal = 1 if boundary == 'top' else -1
            vz_compensation = (target_flux / abs(self.Lx)) * normal
            
            ranges = [self.x.min(), self.x.max()]
            
            self.scenario._print_verbose(f"Applying compensating Vz = {vz_compensation:.3e} to {boundary}")
            self.set_region(boundary, ranges, vz=vz_compensation, mode='add')
        

        elif boundary in ['left', 'right']:
            normal = 1 if boundary == 'right' else -1
            vx_compensation = (target_flux / abs(self.Lz)) * normal
            
            ranges = [self.z.min(), self.z.max()]
            
            self.scenario._print_verbose(f"Applying compensating Vx = {vx_compensation:.3e} to {boundary}")
            self.set_region(boundary, ranges, vx=vx_compensation, mode='add')
        
        return self

    def _evaluate_bound(self,boundary):
        '''
        Internal method to evaluate which coordinates (x or z) must to be used
        '''
        if boundary in ['left', 'right']:
            varying_coord = self.z
        elif boundary in ['top', 'bot']:
            varying_coord = self.x
        else:
            raise ValueError(f"Unknown boundary: {boundary}")

        return varying_coord

    def _apply_veloc(self, mask, boundary, vx, vz, mode='set'):
        '''
        Internal method to apply the velocity into a range, using a mask and a mode.
        '''
        if mode not in ['set', 'add']:
            raise ValueError(f"Invalid mode '{mode}'. Use 'set' or 'add'.")

        if vx is not None:
            
            target_vx = self.velocs[boundary]['vx'] + vx if mode == 'add' else vx
            
            self.velocs[boundary]['vx'] = xr.where(
                mask, 
                target_vx, 
                self.velocs[boundary]['vx']
            )
            
        if vz is not None:
            target_vz = self.velocs[boundary]['vz'] + vz if mode == 'add' else vz
            
            self.velocs[boundary]['vz'] = xr.where(
                mask, 
                target_vz, 
                self.velocs[boundary]['vz']
            )

        return self

    def view_veloc(self):
        '''
        Function to plot the velocities in the boundaries.
        The first figure represents the left and right, and the second the top and bot boundaries. 

        returns: (plt.Figure, plt.Figure)
        '''
        fig1, axslr = plt.subplots(1,2, sharex=True, sharey=True)

        self.velocs.right.vx.plot(y='z',ax=axslr[1],label='vx')
        self.velocs.right.vz.plot(y='z',ax=axslr[1],label='vz')

        self.velocs.left.vx.plot(y='z',ax=axslr[0],label='vx')
        self.velocs.left.vz.plot(y='z',ax=axslr[0],label='vz')
        axslr[0].set_title('left')
        axslr[1].set_title('right')
        axslr[0].legend()
        axslr[0].set_xlabel('vel.')
        axslr[1].set_xlabel('vel.')

        axslr[0].grid()
        axslr[1].grid()

        fig2, axstb = plt.subplots(2,1, sharex=True, sharey=True)

        self.velocs.top.vx.plot(x='x',ax=axstb[0],label='vx')
        self.velocs.top.vz.plot(x='x',ax=axstb[0],label='vz')

        self.velocs.bot.vx.plot(x='x',ax=axstb[1],label='vx')
        self.velocs.bot.vz.plot(x='x',ax=axstb[1],label='vz')

        axstb[0].set_title('top')
        axstb[1].set_title('bot')
        axstb[0].legend()
        axstb[1].set_ylabel('vel.')
        axstb[0].set_ylabel('vel.')

        axstb[0].grid()
        axstb[1].grid()

        print(f'Int V•dS = {self.integrate_normal_components().values}')

        return fig1,fig2


class TemperatureFieldBuilder:
    
    '''
    Constructs the 2D thermal field for Mandyoc scenarios, using the created spatial grid 
    and material data structures.

    Parameters
    ----------
    scenarioBuilder: MandyocBuilder
        The primary scenario instance containing grid geometry and pre-populated material 
        property fields.

    c_cap: float, optional
        Specific heat capacity (in J/(kg K)). Defaults to 1250.0.

    g: float, optional
        Gravitational acceleration (in m2/s). 
        Default is -10.0.

    alpha: float, optional
        Volumetric thermal expansion coefficient (in K^-1). 
        Default is 3.28e-5.

    Attributes
    ----------

    scenario: MandyocBuilder
        Reference pointer to the parent scenario builder.

    Nx, Nz: int
        Number of computational nodes in the x (horizontal) and z (vertical) directions.

    Lx, Lz: float
        Total physical dimensions of the domain along the x and z axes.

    x, z: array-like
        1D spatial coordinate arrays.

    thick_air: float
        Vertical thickness of the upper sticky air boundary layer.

    g, alpha, c_cap: float 
        Stored physical constants governing thermomechanical behavior.

    rho: xarray.DataArray 
        2D spatial density field mapped from the scenario tree.

    kappa: xarray.DataArray
        2D spatial thermal diffusivity field mapped from the scenario tree.

    H: xarray.DataArray 
        2D spatial radiogenic heat production field mapped from the scenario tree.

    k_cond: xarray.DataArray 
        2D thermal conductivity field, automatically computed upon initialization 
        as the product of kappa, rho, and c_cap.

    temperature: xarray.DataArray 
        2D thermal field array of shape (Nx, Nz) initialized to 0.0. Coordinates are 
        mapped to x and z, with internal attributes explicitly defining units as 'deg C'.

    z_corr: xarray.DataArray
        Vertical coordinates vertically offset by the sticky air thickness (z + thick_air) 
        to evaluate true structural depth from the geological surface.

    XX, ZZ (numpy.ndarray)
        2D spatial coordinate meshgrids generated from x and z_corr. Used for fast 
        vectorized matrix masking, boolean logic, and depth-dependent numerical conditions.
    '''

    def __init__(self, scenarioBuilder:MandyocBuilder, 
                 c_cap:float=1250.0,g:float=-10.0,alpha:float=3.28e-5):


        self.scenario = scenarioBuilder
        
        self.Nx = scenarioBuilder.Nx
        self.Nz = scenarioBuilder.Nz
        self.Lx = scenarioBuilder.Lx
        self.Lz = scenarioBuilder.Lz

        self.x = scenarioBuilder.x
        self.z = scenarioBuilder.z
        self.thick_air = scenarioBuilder.thick_air

        self.g = g
        self.alpha = alpha
        self.c_cap = c_cap
        
        self.rho = scenarioBuilder.DTree.fields.rho
        self.kappa = scenarioBuilder.DTree.fields.k
        self.H = scenarioBuilder.DTree.fields.H
        self.k_cond = self.kappa*self.rho*self.c_cap
        
        self.temperature = xr.DataArray(np.zeros((self.Nx,self.Nz),np.float64),dims=('x','z'),
                                     coords={'x':self.x,'z':self.z})
        
        self.temperature.attrs['unit']='deg C'
        
        self.z_corr = (self.temperature.z + self.thick_air)
        self.XX, self.ZZ = np.meshgrid(self.x,self.z_corr)


    def apply_basic_temperature(self, lithosphere_thickness:float,t_pot:float=1350.0,):
        '''
        This function initiate the temperature field combining a linear gradient within the 
        lithosphere and an adiabatic gradient in the mantle.

        Future implementations: make it for each column based on the interface names instead 
        assuming a constant lithosphere.

        Parameters
        ----------
        lithosphere_thickness: float
            The vertical thickness of the lithosphere used to calculate the initial 
            linear thermal gradient.

        t_pot: float, optional
            The potential temperature of the asthenosphere/mantle boundary.
            
            Default is 1350.0.
        '''
        self.t_pot = t_pot
        self.lithosphere_thickness = lithosphere_thickness
        self.temperature[:,:] = (t_pot)/lithosphere_thickness * -self.ZZ.T
        temp_adiabatic = t_pot / np.exp(self.g * self.alpha * -self.ZZ.T / self.c_cap)
        

        self.temperature = xr.where(self.temperature < 0, 0.0,self.temperature)
        self.temperature = xr.where(self.temperature>temp_adiabatic,temp_adiabatic,self.temperature)

    def solve_heat_diffusion2D(self, time_max: float, dt_years: float=0.0):
        '''
        This function solves the 2D transient heat diffusion equation over a specified maximum 
        time, using a Numba-optimized finite difference loop.

        This function is very slow for medium-fine grids, it still need improvements.

        Parameters
        ----------
        time_max: float
            The total simulation time in years for the heat diffusion process.

        dt_years: float, optional
            The time step size in years. If it is 0.0 or exceeds the numerical stability 
            (CFL) limit, the function will automatically cap it to a stable maximum.
            
            Default is 0.0.
        '''
        import numpy as np
        
        # Convert dt from years to seconds to match SI units
        
        dt = dt_years * SEC_PER_YEAR
        num_steps = int(time_max / dt_years)

        # Extract raw numpy arrays from xarray DataArrays for computation speed
        T = self.temperature.transpose('x', 'z').values
        kappa = self.kappa.transpose('x', 'z').values
        H = self.H.transpose('x', 'z').values
        
        # Grid spacing
        dx = abs(self.x[1] - self.x[0])
        dz = abs(self.z[1] - self.z[0])
        
        # Ensure numerical stability (CFL condition)
        kappa_max = np.max(kappa)
        dt_max = (dx**2 * dz**2) / (2 * kappa_max * (dx**2 + dz**2))            

        cond = (-self.ZZ.T >= (-self.lithosphere_thickness)) & (-self.ZZ.T <= 0)

        if (dt > dt_max) or (dt_years<=0):
            
            # Automatically cap the timestep to the maximum stable limit to prevent explosion
            dt = dt_max * 0.95
            num_steps = int((time_max * SEC_PER_YEAR) / dt)

        self.scenario._print_verbose(f"using dt = {dt/SEC_PER_YEAR:.2e} yrs (dt_max is {dt_max/SEC_PER_YEAR:.2e}).")
        self.scenario._print_verbose(f"{num_steps} steps to run.")

        T = _numba_diffusion_loop(
            T, kappa, H, self.c_cap, dx, dz, dt, num_steps, cond
        )

        # Push updated values back to the DataArray
        if self.temperature.dims == ('z', 'x'):
            self.temperature.values = T.T
        else:
            self.temperature.values = T


    def solve_heat_diffusion1D(self, time_max: float, dt_years: float=0.0, x_inx:int=0):
        '''
        This function solves the 1D transient heat diffusion equation on a single vertical 
        column and replicates the resulting thermal profile across the entire 2D domain.

        Future implementations: Apply this calculated profile in a determined range instead in 
        the whole scenario. However, it is possible to create more than one temperature builder
        and to filter the thermal field, applying both builders in differente regions.  

        Parameters
        ----------
        time_max: float
            The total simulation time in years for the 1D heat diffusion process.

        dt_years: float, optional
            The time step size in years. If it is 0.0 or exceeds the numerical stability limit, 
            the function will automatically compute and apply a stable time step.
            
            Default is 0.0.

        x_inx: int, optional
            The index along the x-axis from which the 1D vertical column is extracted 
            for the calculation.
            
            Default is 0.
        '''
        
        dt = dt_years * SEC_PER_YEAR
        num_steps = int(time_max / dt_years)
        print(f"Running {num_steps} iterations for 1D profile...")
        
        # 1. Extract 1D profiles (using the first column at x=0)
        # We use .isel() to safely slice the xarray without assuming axis order
        T_1d = self.temperature.isel(x=x_inx).values.copy()
        kappa_1d = self.kappa.isel(x=x_inx).values.copy()
        H_1d = self.H.isel(x=x_inx).values.copy()

        cond = (self.z_corr.values >= (-self.lithosphere_thickness)) & (self.z_corr.values < 0)
        # print(f'z>{(-self.lithosphere_thickness)}')
        # print(f'z<{0}')
        # print(self.z_corr[cond])
        
        dz = abs(self.z[1] - self.z[0])
        
        # 2. CFL Condition (1D limit is less restrictive than 2D)
        kappa_max = np.max(kappa_1d)
        dt_max = (dz**2) / (2 * kappa_max)
        
        if dt > dt_max:
            print(f"Warning: Input dt exceeds 1D stable limit ({dt_max/SEC_PER_YEAR:.2e} yrs).")
            dt = dt_max * 0.99
            num_steps = int((time_max * SEC_PER_YEAR) / dt)

        self.scenario._print_verbose(f"using dt = {dt/SEC_PER_YEAR:.2e} yrs (dt_max is {dt_max/SEC_PER_YEAR:.2e}).")
        self.scenario._print_verbose(f"{num_steps} steps to run.")
        
        # 3. 1D Finite Difference Loop
        for step in range(num_steps):
            T_new = np.copy(T_1d)
            
            # Derivatives along the Z-axis
            dT_dz = (T_1d[2:] - T_1d[:-2]) / (2 * dz)
            dK_dz = (kappa_1d[2:] - kappa_1d[:-2]) / (2 * dz)
            d2T_dz2 = (T_1d[2:] - 2 * T_1d[1:-1] + T_1d[:-2]) / (dz**2)
            
            # Diffusion term
            diffusion_z = kappa_1d[1:-1] * d2T_dz2 + dK_dz * dT_dz
            
            # Forward Euler update
            T_new[1:-1] = T_1d[1:-1] + dt * (diffusion_z + H_1d[1:-1] / self.c_cap)
            
            # Boundary Conditions (Dirichlet: Fixed temperatures at top and bottom)
            T_new[0] = T_1d[0]
            T_new[-1] = T_1d[-1]
            
            T_1d = np.where(cond, T_new, T_1d)
            
        # 4. Replicate and map back to 2D
        nx = len(self.x)
        
        # Verify the dimensions of the host DataArray to broadcast correctly
        if self.temperature.dims == ('z', 'x'):
            # Tile T_1d into columns: shape becomes (nz, nx)
            T_2d = np.tile(T_1d[:, np.newaxis], (1, nx))
        else:
            # Tile T_1d into rows: shape becomes (nx, nz)
            T_2d = np.tile(T_1d, (nx, 1))
            
        self.temperature.values = T_2d