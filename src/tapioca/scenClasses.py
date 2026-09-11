import os, gc, json, glob
from pathlib import Path
from shapely.geometry import Polygon
from shapely import to_ragged_array

from scipy.interpolate import RegularGridInterpolator

import numpy as np
import xarray as xr
from xarray import DataTree

from ._variables import VARS_TYPES, VARIABLES_LIST, INTERFACES_PARAMETERS,DEFAULT_MATERIAL,MATERIAL_PARAMETERS
from ._aux_functions import read_params

#Mandyoc Scenario class
class MandyocScen:
    """
    Core data manager for Mandyoc simulations (results).

    This class parses output directories, reads model parameters, and constructs 
    a unified DataTree. The classes alsos contains the params.txt stored as a
    frozen dictionary. 
    
    It aligns meshes (common variables), upscaled lithology, surface topography,
    and Lagrangian particles into a single navigable object. 

    Functions to manipulate data of `MandyocScen` are built-in within this class
    and as acessors in `post_processing.py`.

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
    set_nan : float, optional
        Value to be considered "NoData" in X and Z. The values will be 
        converted to nan. If None, `set_nan` is not considered.

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
    XMAX,XMIN,ZMIN,ZMAX : float
        The domain boundaries for the loaded scenario. They are equivalent to the Lx 
        and Lz parameters if the whole domain is loaded.
    dx,dz : float
        The node's spacing in x and z (constants).
    dt : np.ndarray
        The time spacing between the steps, which varies over the simulation.
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
                 filter_air:bool=True, air_layer:int=None, set_nan:float=None, #only relevant if load particles
                 verbose:bool=False):

        # Setting directories and scen name
        self.path = Path(path)
        self.params = read_params(os.path.join(path,'param.txt'))
        self.verbose = verbose
        if isinstance(name, str): self.name = name
        else: self.name = self.path.name

        self._print_verbose(f'Scenario at: {self.path}')
        self._print_verbose(f'Scenario name: {self.name}')
        self._print_verbose(f'Params - Ok')
        
        # Handle Variables and Extract Metadata
        if isinstance(variables, str): variables = [variables]
        elif not isinstance(variables, (list, tuple, np.ndarray)): variables = list(variables)
        
        if len(variables) == 0:
            raise ValueError("The 'variables' list cannot be empty. Need at least one variable to extract metadata.")

        self.dx = None
        self.dz = None
        self.dt = None
        self.get_scenarioData(variables[0])
        
        self.xlimits = xlimits if xlimits is not None else [self.XMIN, self.XMAX]
        self.zlimits = zlimits if zlimits is not None else [self.ZMIN, self.ZMAX]
        self.tlimits = tlimits if tlimits is not None else [self.TMIN, self.TMAX]
        self.thick_air = thick_air #m

        
        self._print_verbose(f'x limits: {self.xlimits} (dx = {self.dx: 11.4f})')
        self._print_verbose(f'z limits: {self.zlimits} (dz = {self.dz: 11.4f})')
        self._print_verbose(f'time limits: {self.tlimits} (average dt = {np.mean(self.dt): 3.8f})')
        
        self.z_corrected = False
        self.particles_loaded = False
        
        # Initialize the empty DataTree
        self.DTree = DataTree.from_dict( {"mesh":None,"surface":None,"particles":None})

        # Passing some metadata to the DataTree
        # self.DTree.attrs['name'] = self.name
        # self.DTree.attrs['xlimits'] = self.xlimits
        # self.DTree.attrs['zlimits'] = self.zlimits
        # self.DTree.attrs['tlimits'] = self.tlimits
        
        # Reading and storing standard variables (Nx,Nz)
        standard_datasets = []
        for var in variables:
            ds = self._load_spatial_var(var, chunks=chunks_vars)
            if ds is not None:
                standard_datasets.append(ds)
        
        if standard_datasets:
            self.DTree['/mesh/original'] = xr.merge(standard_datasets, join='outer')

        self._print_verbose(f"Variables loaded: {' '.join(variables)}")
        
        # Reading and storing the lithology (upscaled mesh)
        if load_lithology==True:
            ds_litho = self._load_spatial_var('lithology', chunks=chunks_vars)
            if ds_litho is not None:
                self.DTree['/mesh/upscaled'] = ds_litho

        self._print_verbose("Lithology loaded")
            
        # Reading the surface/topography (only X dimension)
        if load_surface==True:
            
            ds_surf = self._load_spatial_var('surface', chunks=chunks_vars)
            if ds_surf is not None:
                self.DTree['/surface/topography'] = ds_surf

        self._print_verbose("Surface loaded")
            
        # Loading particles
        if load_particles:
            self._load_particles(particles_file, chunks={'id': 'auto'}, filter_air=filter_air, air_layer=air_layer,set_nan=set_nan)
            
        self._print_verbose(f"Particles [{particles_file}] loaded")
        
        return None

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

            self.dx = abs(float(ds.x[1] - ds.x[0]))
            self.dz = abs(float(ds.z[1] - ds.z[0]))
            self.dt = np.diff(ds.time.values)

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
            self._print_verbose("Z was already corrected")
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
        
        
        self._print_verbose(f"Z coordinate corrected by + {factor} m")
        self._print_verbose(f"New z limits: {self.zlimits}")
            
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
        if 'step' in v.data_vars:
            v = v.set_coords('step')
        
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
    
    def _load_particles(self, name:str, chunks:dict={'id': 'auto'}, filter_air:bool=True, 
                            air_layer:int=None, set_nan:float=None):
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
        
        set_nan : float, optional
            Value to be considered "NoData" in X and Z (e.g. -999.). The values will be 
            converted to nan. If None, `set_nan` is not considered.

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
        if "step" in particles.data_vars:
            particles = particles.set_coords('step')
        
        particles = particles.sel(time=slice(self.tlimits[0], self.tlimits[-1]))

        if set_nan is not None:
            particles['x'] = particles['x'].where(particles['x'] != set_nan)
            particles['z'] = particles['z'].where(particles['z'] != set_nan)

        if filter_air and 'lithology' in particles.data_vars:
            if isinstance(air_layer, int): air = air_layer
            else: air = int(particles.lithology.max()) #future: find air layer by density?
                
            cond = particles.lithology != air
            particles = particles.where(cond)
            self._print_verbose(f'Air particles filtered [{air}]')
            
        self.DTree['/particles/original'] = particles
        self.particles_loaded = True
        return True

    
    def _apply_selection(self, valid_ids:np.ndarray, 
                         source_name:str='', target_name:str=''):
        """
        Internal function to select a dataset by the IDs and store it in the DataTree.
        """

        paths, names, pts = self._eval_particles_selection(source_name=source_name,target_name=target_name,get_selection=True)
        source_path, target_path = paths
        source_name, target_name = names

        self._print_verbose(f'Selected particles: {source_name}')
        self._print_verbose(f'Selection at: {target_name}')

        # Lazy selection
        pts = self.DTree[source_path].ds
        self.DTree[target_path] = pts.sel(id=valid_ids)
    
        gc.collect()
        return None

    def _eval_particles_selection(self,source_name:str, target_name:str, 
                                  get_selection:bool=False):
        """
        Internal method that evaluates and validates particle node paths and names for selection operations.

        Parameters
        ----------
        source_name : str
            The input node name.
        target_name : str
            The desired output node name.
        get_selection : bool, optional
            If True, returns the dataset object alongside paths and names.

        Returns
        -------
        tuple
            (paths, names, dataset) or (paths, names) where paths and names are tuples 
            of (source, target).
        """

        if len(target_name) == 0:
            target_name = 'selected'
        
        elif len(source_name) == 0:
            source_name = 'original'
        
        source_path = f'/particles/{source_name}'
        target_path = f'/particles/{target_name}'

        if get_selection:
            return (source_path,target_path), (source_name,target_name), self.DTree[source_path].ds
        else:
            return (source_path,target_path), (source_name,target_name)
    
    #Selecting particles
    def selectParticles_bytimerange(self, timerange:list, 
                    source_name:str='', target_name:str=''):
        """
        Selects particles that entered the domain (i.e., appeared) between two timesteps.

        This identifies particles present at `timerange[1]` that were not present 
        at `timerange[0]` to create a selection.

        Parameters
        ----------
        timerange : list of float
            Time range [t_start, t_end] in Myr.
        source_name : str, optional
            Node name in `'/particles/'` to filter from.
        target_name : str, optional
            Node name in `'/particles/'` where the selection will be saved.

        Returns
        -------
        MandyocScen
            Returns the class instance for method chaining.
        """
        # apply the support selection function
        paths, names,pts = self._eval_particles_selection(source_name=source_name,target_name=target_name,get_selection=True)
        source_path, target_path = paths
        source_name, target_name = names
        
        # select particles based on two snapshots
        snap_0 = pts.sel(time=timerange[0], method='nearest').compute()
        snap_i = pts.sel(time=timerange[1], method='nearest').compute()

        # selecting particles based on the X coordinate and discarding null values
        tr_0 = snap_0.id.where(~snap_0.x.isnull(), drop=True).values
        tr_i = snap_i.id.where(~snap_i.x.isnull(), drop=True).values

        # get the difference between them (using sets)
        ids = list(set(tr_i)-set(tr_0))

        self._print_verbose(f'{len(ids)} particles were selected between {timerange[0]}-{timerange[1]}')
            
        # apply the support selection function
        self._apply_selection(ids, 
                              source_name=source_name,
                              target_name=target_name)
        return self
        
    
    def selectParticles_bycoords(self, xlim:list=None, zlim:list=None, tsel:float=None, 
                source_name:str='',target_name:str=''):
        """
        Selects particles within a specific spatial box at a given timestep.

        Parameters
        ----------
        xlim : list of float, optional
            Horizontal bounds [xmin, xmax]. If None, uses `self.xlimits`.
        zlim : list of float, optional
            Vertical bounds [zmin, zmax]. If None, uses `self.zlimits`.
        tsel : float, optional
            The specific time (Myr) to perform the spatial selection. 
            If None, defaults to the first timestep (0).
        source_name : str, optional
            Node name in `'/particles/'` to filter from. Defaults to 'original'.
        target_name : str, optional
            Node name in `'/particles/'` where the selection will be saved.

        Returns
        -------
        MandyocScen
            Returns the class instance for method chaining.
        """
        
        if xlim is None: xlim = self.xlimits
        if zlim is None: zlim = self.zlimits

        paths, names, pts = self._eval_particles_selection(source_name=source_name,target_name=target_name,get_selection=True)
        source_path, target_path = paths
        source_name, target_name = names
            
        if tsel is None: tsel = 0
        
        pts = pts.sel(time=tsel, method='nearest')
        condX = ((pts["x"] >= xlim[0]) & (pts["x"] <= xlim[1])).compute()
        pts =  pts.where(condX, drop=True)
        
        condZ = ((pts["z"] >= zlim[0]) & (pts["z"] <= zlim[1])).compute()
        pts =  pts.where(condZ, drop=True)
        ids = pts.id.values
              
        # apply the support selection function
        self._apply_selection(ids, 
                              source_name=source_name,
                              target_name=target_name)
        
        return self
    
    
    def selectParticles_bylithology(self, lithologies:list, tsel:float=None, source_name:str='',target_name:str=''):
        """
        Selects particles based on their lithological phase at a given timestep.

        Parameters
        ----------
        lithologies : list of int
            List of lithology IDs to select.
        tsel : float, optional
            The specific time (Myr) to perform the selection. 
            If None, defaults to 0.
        source_name : str, optional
            Node name in `'/particles/'` to filter from. Defaults to 'original'.
        target_name : str, optional
            Node name in `'/particles/'` where the selection will be saved.

        Returns
        -------
        MandyocScen
            Returns the class instance for method chaining.
        """
        
        paths, names,pts0 = self._eval_particles_selection(source_name=source_name,target_name=target_name,get_selection=True)
        source_name, target_name = names
        pts = pts0.copy()

        if tsel is None: tsel = 0 # Future: to create a function to evaluate an automatic tsel
        
        pts = pts.sel(time=tsel, method='nearest')
        cond = pts.lithology.isin(lithologies).compute()
        ids =  pts.id.where(cond, drop=True).values
        
        self._apply_selection(ids, 
                              source_name=source_name, 
                              target_name=target_name)
        return self
    
    
    def classify_ParticlesRange(self, domain_intervals:dict, tsel:float=None,
                               source_name:str='',target_name:str=''):
        
        """
        Classifies all particles into domains based on X-coordinate ranges at a specific timestep.

        The classification result is stored as a new coordinate variable in the particles dataset.

        Parameters
        ----------
        domain_intervals : dict
            Dictionary defining categories and their X-ranges. Must include 
            a 'field_name' key (e.g., {'field_name': 'my_domain', 'necking': [100e3,200e3], 'distal': [200e3, 300e3]}).
        tsel : float, optional
            The specific time (Myr) for the classification snapshot.
        source_name : str, optional
            Node name in `'/particles/'` to classify.
        target_name : str, optional
            Node name in `'/particles/'` where the updated dataset is saved.

        Returns
        -------
        MandyocScen
            Returns the class instance for method chaining.
        """
        #Classify all particles based on X ranges, given a time step tsel
        #Categories are based on the domain intervals keys
        
        paths, names,pts0 = self._eval_particles_selection(source_name=source_name,target_name=target_name,get_selection=True)
        source_path, target_path = paths
        source_name, target_name = names
        pts = pts0.copy()
        
        if tsel is None: tsel = 0
        
        domain_intervals = domain_intervals.copy()
        
        try: field_name = domain_intervals['field_name']
        except: field_name = 'domain'
        del domain_intervals['field_name']
        
        typename = type(list(domain_intervals.keys())[0])
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
        
        target_path = f'/particles/{target_name}'
        self.DTree[target_path] = pts
        
        pts[field_name].attrs['reference timestep'] = f'{tsel}myr'
        pts[field_name].attrs['classes range'] = str(domain_intervals)
        gc.collect()
        
        return self

    def fieldToParticle(self, variable:str, source_name:str='', method:str='linear'):
        """
        Interpolates an Eulerian mesh field onto Lagrangian particle trajectories.

        This method extracts a specified spatial field from the original mesh 
        node (`self.DTree.mesh.original`) and interpolates its values at the 
        spatiotemporal coordinates (`x`, `z`, `time`) of the given particles. 
        The resulting interpolated variable is then appended directly to the 
        particle dataset.

        Parameters
        ----------
        variable : str
            The name of the variable from the Eulerian mesh to interpolate 
            (e.g., 'temperature', 'pressure').
        source_name : str, optional
            The name of the selected particles to update (e.g., 
            'original', 'selected', etc). If left empty, the internal evaluation 
            will select 'original'.
            
            Default is `''`.
        method : str, optional
            The interpolation method to be used by `xarray.Dataset.interp`. 
            Common options include 'linear' and 'nearest'.
            
            Default is `'linear'`.

        Returns
        -------
        MandyocScen
            Returns the class instance (`self`) to allow for method chaining.

        Notes
        -----
        **Tree Mutation:** This method modifies the `DataTree` in-place by 
        assigning the newly interpolated variable to the dataset located at 
        the evaluated `source_path` node.
        
        **Memory Management:** The method explicitly deletes temporary datasets 
        and triggers garbage collection (`gc.collect()`) to prevent memory spikes 
        during interpolation.
        """
        paths, _, pts0 = self._eval_particles_selection(source_name=source_name,target_name='',get_selection=True)
        source_path, _ = paths

        pts = pts0.copy()

        field = self.DTree.mesh.original[variable]
        
        field_interpolated = field.interp(
            x=pts.x, 
            z=pts.z,
            time=pts.time,
            method=method
        )
    
        self.DTree[source_path]=self.DTree[source_path].assign({variable:field_interpolated})

        del pts, field_interpolated, pts0
        gc.collect()

        return self


class MandyocBuilder:

    def __init__(self, path:str, Nx:int, Nz:int, Lx:float, Lz:float, Ny:int=0, Ly:float=0, verbose:bool=True):
        self.path = path
        self.verbose = verbose

        self.Nx = Nx
        self.Nz = Nz
        self.Lx = Lx
        self.Lz = -Lz
        self.Ly = self.Ny = self.y = None
        self.dimensions = 2

        self.x = np.linspace(0, self.Lx, self.Nx)
        self.z = np.linspace(self.Lz, 0, self.Nz)

        # Creating the DataTree structure to store the model layers
        self.DTree = DataTree.from_dict( {"field":
                            xr.Dataset(coords={'x': self.x, 'z': self.z}), 
                                    "interfaces":xr.Dataset(coords={'x': self.x})})

        self.DTree['/field']["x"].attrs["units"] = "m"
        self.DTree['/field']["z"].attrs["units"] = "m"
        self.DTree['/field'].attrs["description"] = "This group contains the field of the model, each variable is a DataArray considering the material properties from the interfaces."

        self.DTree['/interfaces']["x"].attrs["units"] = "m"
        self.DTree['/interfaces'].attrs["description"] = "This group contains the interfaces of the model, each interface is a DataArray with the depth and the material properties as attributes."
        
        # setting the model dimensions
        if Ny>0 and abs(Ly)>0:
            self.Ly = Ly
            self.Ny = Ny
            self.dimensions = 3
            self.y = np.linspace(0, self.Ly, self.Ny)
            self.DTree = DataTree.from_dict( {"field":
                                xr.Dataset(coords={'x': self.x, 'y': self.y, 'z': self.z}), 
                                        "interfaces":xr.Dataset(coords={'x': self.x, 'y': self.y})})

            self.DTree['/field'].expand_dims(dim='y', axis=1).assign_coords(y=self.y)
            self.DTree['/interfaces'].expand_dims(dim='y', axis=1).assign_coords(y=self.y)
            
            self.DTree['/field']["y"].attrs["units"] = "m"
            self.DTree['/interfaces']["y"].attrs["units"] = "m"


    def create_interface(self, id_layer:int, position, interface_name:str=None, **material):

        if interface_name==None:
            interface_name = f'interface_{id_layer}'
                
        # Defining material properties
        interface_parameters = DEFAULT_MATERIAL.copy() # default material

        for k in list(material.keys()):
            interface_parameters[k] = material[k]

            if self.verbose==True:
                print(f'{k}:{interface_parameters[k]}')

        if interface_parameters['rho'] == None:
            raise ValueError("ERROR: You must to set a density (rho) to the material interface!")


        if id_layer!=-1:
            interface_coords = self._evaluate_interface_position(position)
            #create a default interface at the top of the model (ZMAX)

        else:
            print(f'''Layer with ID {id_layer} is the top layer, it contains the initial material 
            properties of the model that is above the last interface.\n
            You can set material properties for this layer, but you cannot create an interface.''')
            
            interface_coords = np.ones(self.Nx) if self.dimensions==2 else np.ones((self.Nx,self.Ny))
            interface_coords *= 0

        if self.dimensions==2:
            self.DTree['/interfaces'][interface_name] = (xr.DataArray(interface_coords, coords=[self.x], 
                                                                              dims=['x'], name=interface_name))
            

        elif self.dimensions==3:
            self.DTree['/interfaces'][interface_name] = (xr.DataArray(interface_coords, coords=[self.x,self.y], 
                                                                                 dims=['x','y'], name=interface_name))


        self.DTree['/interfaces'][interface_name].attrs['id'] = id_layer
        self.DTree['/interfaces'][interface_name].attrs.update(interface_parameters)
        return self


    def add_interface_from_shapely(self, layer:str, polygon:Polygon, poly_name:str=None, **material): # 2D
    
        polycoords = to_ragged_array([polygon])[1]
        xmin,zmin, xmax,zmax = polygon.bounds
        mid_z = (zmin+zmax)/2

        if self.verbose:
            print(f'Polygon bounds: ({xmin,zmin}) ({xmax,zmax})')
            print(f'Mid line is at {mid_z} m')
        
        top_args = polycoords[:,1] >= mid_z
        bot_args = polycoords[:,1] < mid_z
        
        x_top = np.append(0, polycoords[:,0][top_args][1:])
        x_top = np.append(x_top, self.Lx)

        z_top = np.append(mid_z,polycoords[:,1][top_args][1:])
        z_top = np.append(z_top, mid_z)

        x_bot = np.append(0, polycoords[:,0][bot_args][1:])
        x_bot = np.append(x_bot, self.Lx)

        z_bot = np.append(mid_z,polycoords[:,1][bot_args][1:])
        z_bot = np.append(z_bot, mid_z)

        bot_interface = np.array([x_bot,z_bot]).T
        top_interface = np.array([x_top,z_top]).T

        print(bot_interface)
        
        # handle with the ids
        interfaces_names = list(self.DTree.interfaces.variables)[:-1]
        id_above = self.DTree['/interfaces'][layer].attrs['id']

        for n in interfaces_names:
            if self.DTree['/interfaces'][n].attrs['id'] == id_above:
                above_material = self.DTree['/interfaces'][n].attrs.copy()
                if self.verbose:
                    print(f"Above layer is: ({self.DTree['/interfaces'][n].attrs['id']}) {n}")

            if self.DTree['/interfaces'][n].attrs['id'] >= id_above:
                self.DTree['/interfaces'][n].attrs['id'] += 2


        self.create_interface(id_above,bot_interface,interface_name=f'bot_{poly_name}',**above_material)
        self.create_interface(id_above+1,top_interface,interface_name=f'top_{poly_name}',**material)


        return self

    def _evaluate_interface_position(self, coord):
        """
        internal method to verify what type of coordinate or position the interface
        is recieving 
        """

        if isinstance(coord,float) or isinstance(coord,int):
            print('The interface is a flat line/surface')
            coordArray = np.ones(self.Nx) if self.dimensions==2 else np.ones((self.Nx,self.Ny))
            return coordArray*coord

        # Verify if the coord is a bunch of points to interpolate or is the whole interface (line/surface)
        elif isinstance(coord,np.array) and self.dimensions==2:
            print('The interface is a line, interpolating the points to create the interface')

            if len(coord[0])!=2:
                raise ValueError("ERROR: The coordinates for the interface must be a 2D array with two columns (x,z).")

            args_sorted = np.argsort(coord[:,0])
            
            return np.interp(self.x, coord[:,0][args_sorted], coord[:,1][args_sorted])
        
        else:
            raise ValueError("ERROR: The coordinates for the interface must be a float, int for 2D and 3D models, or a 2D array with two columns (x,z) for 2D models.")

        # Add interfaces from Shapely geometries
        
        return None
    
    def set_params(self, **kwargs):
        """
        Set parameters for the model
        """

        # create a function to evaluate params!

        self.params = kwargs
        
        return self