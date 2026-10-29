'''
Builder of the input files of a Mandyoc scenario drawn in an SVG file.

The geometry and the materials come from the SVG (see tapioca.inputs.svg.reader). The shape
labelled "setup" gives the domain, the resolution and, optionally, the numerical
experiment:

    Lx = 1600 km
    Lz = 400 km
    dx = 2 km
    dz = 2 km
    regime = accordion                         (rifting, compression, accordion, quiescence)
    velocity = 1 cm/yr
    moving_boundary = right                    (left, right or symmetric)
    quiescence_time = 0 Myr
    time_max = 500 Myr
    magmatism = off
    magmatism_extraction = off
    magmatism_Ts0 = 1080
    preset_initial_temperature = ../preliminary
    temperature_solver = 2D                    (2D or columns)
    min_thickness = 1 km

Each setting may also be given to SVGScenario (or on the command line, tapioca-svg),
which overrides the SVG; without it, the default is used.

Example
-------
    from tapioca.inputs.svg import SVGScenario
    scenario = SVGScenario("rifting.svg", regime="accordion", quiescence_time=0)
    scenario.write_all("my_scenario")

The scenario is built on the generic builders of tapioca: scenario.builder is a MandyocBuilder with the
interfaces and the material fields of the SVG, and build_temperature and build_velocity use a
TemperatureFieldBuilder and a VelocityFieldBuilder. The input files are exported by scenario.builder,
so interfaces.txt has the 13 parameters of tapioca.DEFAULT_MATERIAL.
'''

import os, glob, shutil, platform

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from ..._variables import SEC_PER_YEAR, PARAMS_TYPES
from ..._aux_functions import find_nearest, read_preset_temperature, adapt_preset_temperature
from ..generic_input_classes import TemperatureFieldBuilder, VelocityFieldBuilder, VerboseExtension
from ..checks import SVGChecks
from ._variables import (REGIMES, MANTLE_RHEOLOGIES, SYSTEMS, DEFAULT_PARAMS, ORGANIZE_SCRIPT, RUN_SCRIPT,
                         RUN_LAUNCHER, ZIP_SCRIPT)
from .materials import MandyocLayer, BASALT
from ... import __version__
from .reader import (SVGInputError, SETUP_LABEL, MIN_THICKNESS, EXPERIMENT_KEYS, SVG_VERSION, INKSCAPE_VERSION,
                     svg_versions, read_setup, layers_from_svg, layers_to_builder, points_from_spacing, parse_quantity,
                     layer_indices, slot_thickness, layer_thickness)

__all__ = ["SVGScenario", "REGIMES"]

###############################################################################
# Scenario builder
###############################################################################

class SVGScenario(SVGChecks, VerboseExtension):
    '''
    Mandyoc scenario drawn in an SVG file.

    The geometry, the materials and the domain come from the SVG (see tapioca.inputs.svg.reader). The settings
    of the numerical experiment are taken from the arguments, then from the description of the setup
    shape of the SVG and finally from the defaults below.

    Parameters
    ----------
    svg_file : str
        SVG file with the scenario.
    dx, dz : float, optional
        Horizontal and vertical resolution [km] (default: from the SVG).

    Settings of the numerical experiment (keyword arguments, see reader.EXPERIMENT_KEYS); each one overrides the
    setting of the setup shape of the SVG. Numbers may be given with units, as in the SVG (e.g. '1 cm/yr'):

    regime : str, optional
        'rifting' (extension), 'compression' (convergence), 'accordion' (extension, convergence,
        quiescence and extension again, as in M300) or 'quiescence' (no velocity at the boundaries).
        Default: rifting.
    velocity : float, optional
        Velocity at the moving boundary [cm/yr] (with symmetric, the sum of the velocities of both boundaries).
        Default: 1 cm/yr.
    moving_boundary : str, optional
        'left', 'right' or 'symmetric' (rifting and accordion: both boundaries move with velocity/2 and the
        reference frame is at the centre of the domain). Default: right.
    quiescence_time : float, optional
        Duration of the quiescence [Myr] in the accordion and quiescence regimes. Default: 0 Myr.
    time_max : float, optional
        Duration of the simulation [Myr]. Default: 60 Myr in rifting and compression, the end of the
        second extension in accordion (which time_max sets when given) and quiescence_time (or 1000 Myr)
        in quiescence.
    magmatism, magmatism_extraction : str, optional
        'on' or 'off'. Default: off. With both on, a basalt layer of zero thickness receives the extracted
        melt (see magmatic_layer).
    magmatism_Ts0 : float, optional
        magmatism_Ts0 of param.txt [oC]. Default: 1080.
    water_layer : str or int, optional
        Label (or ID from the bottom) of the unit drawn as water, right below the sticky air. The surface (top of
        the solid layers) is then its base. Default: none (no water).
    magmatic_layer : str or int, optional
        Label (or ID from the bottom) of the layer that receives the melt extracted by Mandyoc (magmatic_layer of
        param.txt). Default: with magmatism and magmatism_extraction on, a basalt layer of zero thickness added at
        the top of the solid layers (below the water if there is water, otherwise below the sticky air); otherwise
        the top solid layer. A layer given here is used as it is, so no basalt layer is added.
    preset_initial_temperature : str, optional
        Folder (with temperature.nc, or param.txt and temperature_*.txt) or NetCDF file of a preliminary
        scenario (see tapioca.read_preset_temperature): the horizontal mean of its final
        temperature, adapted to the lithosphere of the SVG, is the initial temperature. A relative path in the
        SVG is relative to the folder of the SVG. Default: none (conductive model).
    temperature_solver : str, optional
        Steady-state conduction in the lithosphere of the conductive model: '2D' (the whole lithosphere at once,
        so heat also flows sideways where the layers change along x) or 'columns' (each column alone).
        Default: 2D.

    Other parameters:

    asthenosphere_label, seed_label : str
        Labels of the asthenosphere and of the weak seed in the SVG.
    sampling : float
        Vertical sampling used to locate the interfaces in the SVG [m].
    min_thickness : float, optional
        Pieces of a unit thinner than min_thickness [km] in every column, when the unit is thicker
        elsewhere, are removed (slivers left by overlapping shapes; 0 keeps them). Default: from the setup
        shape of the SVG, or 1 km.
    verbose : bool
        Print the steps.

    Attributes
    ----------
    layers : list of MandyocLayer (bottom to top); layer.interface is the depth of the layer base [m]
    x, z : coordinates of the nodes [m] (z is the depth, from the bottom to the top of the domain)
    T : initial temperature (Nz, Nx), after build_temperature()
    VX, VZ : initial velocity (Nz, Nx), after build_velocity()
    builder : MandyocBuilder (tapioca.inputs) with the interfaces and the fields of the scenario; the
              input files are exported by it, so it can be changed by script before write_all
    temperature_builder, velocity_builder : TemperatureFieldBuilder and VelocityFieldBuilder of the scenario,
              after build_temperature() and build_velocity()
    '''

    #Initial temperature (the thermal diffusivity k of each layer comes from its description or rheology)
    TP = 1350.0 #mantle potential temperature [oC], at the LAB
    ccapacity = 1250.0 #heat capacity [J/kg/K]

    #Accordion regime [Myr]
    ti_convergence = 25.0
    ti_quiescence = 75.0
    dt_rifting2 = 60.0

    #Velocity boundary conditions
    fac_air = 10.0e3 #no velocity in the sticky air above thickness_sa - fac_air [m]
    h_v_reference = 'lithosphere' #or 'lithospheric_mantle': reference thickness of the region with constant velocity

    #Plots: profiles of temperature and YSE
    min_profile_width = 0.05 #fraction of Lx
    profile_colors = ['tab:orange', 'tab:purple', 'tab:green', 'tab:blue', 'r', 'k'] #the thickest lithosphere is the last one

    def __init__(self, svg_file, dx=None, dz=None, asthenosphere_label='asthenosphere', seed_label='weak_seed',
                 sampling=100.0, min_thickness=None, verbose=True, **settings):

        self.verbose = verbose
        self.svg_file = self._check_file('SVG file', svg_file)
        self._check_keys('setting', settings, EXPERIMENT_KEYS)
        setup = self._check_setup(read_setup(svg_file), svg_file, SETUP_LABEL)
        self.svg_version, self.inkscape_version = svg_versions(svg_file)
        self._check_svg_versions(svg_file, self.svg_version, self.inkscape_version, SVG_VERSION, INKSCAPE_VERSION)

        #domain and resolution
        self.Lx, self.Lz = setup['Lx'], setup['Lz']
        self.Nx = points_from_spacing(self.Lx, dx * 1.0e3, 'dx') if(dx) else setup['Nx']
        self.Nz = points_from_spacing(self.Lz, dz * 1.0e3, 'dz') if(dz) else setup['Nz']
        self._check_domain(self.Lx, self.Lz, self.Nx, self.Nz)
        self.x = np.linspace(0, self.Lx, self.Nx)
        self.z = np.linspace(self.Lz, 0, self.Nz)
        self.X, self.Z = np.meshgrid(self.x, self.z)
        self.dx = self.Lx / (self.Nx - 1)
        self.dz = self.Lz / (self.Nz - 1)

        #numerical experiment: argument > SVG > default
        self._experiment = setup.get('experiment', {})
        self.source = {}
        self._read_settings(settings)

        #layers
        self.asthenosphere_label = asthenosphere_label
        self.seed_label = seed_label
        self.min_thickness = min_thickness * 1.0e3 if(min_thickness is not None) else setup.get('min_thickness', MIN_THICKNESS)
        self.layers = layers_from_svg(svg_file, self.Lx, self.Lz, self.Nx, dz=sampling, setup_id=setup['id'],
                                      min_thickness=self.min_thickness, verbose=verbose)
        self._check_layers(self.layers, svg_file, asthenosphere_label)

        #water and magmatic layers: chosen by the user (label, or ID of the layer from the bottom), never guessed
        self.water_layer = self._layer_index('water_layer', self._setting('water_layer', settings.get('water_layer'), None))
        if(self.water_layer is not None and self.water_layer != len(self.layers) - 2):
            self._fail(f"water_layer '{self.layers[self.water_layer].layer_label}' must be the unit right below the sticky air.")
        magmatic_layer = self._layer_index('magmatic_layer', self._setting('magmatic_layer', settings.get('magmatic_layer'), None))
        if(magmatic_layer is None and self.magmatism == 'on' and self.magmatism_extraction == 'on'):
            #Basalt layer with zero thickness at the top of the solid layers (below the water, if there is water,
            #otherwise below the sticky air): it can not be drawn in the SVG and receives the melt extracted by Mandyoc
            magmatic_layer = self._insert_basalt(len(self.layers) - 1 if(self.water_layer is None) else self.water_layer)
            if(self.water_layer is not None):
                self.water_layer += 1
        #magmatic_layer of param.txt: by default the top of the solid layers (below the water or the sticky air)
        top_solid = len(self.layers) - 2 if(self.water_layer is None) else self.water_layer - 1
        self.magmatic_layer = magmatic_layer if(magmatic_layer is not None) else top_solid

        #sticky air and lithosphere along x
        #surface: top of the solid layers (base of the sticky air, or of the water); its interface ID is the
        #ID of the layer below it (interface i is the top of layer i)
        self.surface_id = len(self.layers) - 2 if(self.water_layer is None) else self.water_layer - 1
        self.thickness_sa = self.layers[self.surface_id + 1].interface #depth of the surface
        self.idx_ast = max(layer_indices(self.layers, asthenosphere_label))
        self.thickness_litho = self.layers[self.idx_ast+1].interface - self.thickness_sa #surface to the top of the asthenosphere
        self.thickness_mantle = np.zeros(self.Nx) #lithospheric mantle (olivine rheologies above the asthenosphere)
        for i in range(self.idx_ast + 1, len(self.layers)):
            if(self.layers[i].rheology_name in MANTLE_RHEOLOGIES):
                self.thickness_mantle += slot_thickness(self.layers, i, slice(None))

        #position of the vertical profile of the scenario infos and the profiles of the plots
        self.x_profile = self.Lx / 4.0
        self.idx_profile = find_nearest(self.x, self.x_profile)
        self.profiles = self._find_profiles()

        #generic builder of tapioca (tapioca.inputs.MandyocBuilder) with the interfaces and the material fields of
        #the SVG; the temperature, the velocity and the YSE are built on it and the input files are exported by it
        self.builder = layers_to_builder(self.layers, self.Lx, self.Lz, Nz=self.Nz, seed_label=seed_label, verbose=False)
        self.builder.create_materials_fields()
        self.temperature_builder, self.velocity_builder = None, None

        self.T, self.VX, self.VZ = None, None, None

    ###########################################################################
    # Settings
    ###########################################################################

    def _setting(self, name, value, default):
        '''value (argument) > SVG > default'''
        if(value is not None):
            self.source[name] = 'argument'
            return value
        if(name in self._experiment):
            self.source[name] = 'SVG'
            return self._experiment[name]
        self.source[name] = 'default'
        return default

    def _layer_index(self, name, value):
        '''Index (ID from the bottom) of the layer given by its label (the top slot of the unit) or its ID; None if value is None.'''
        if(value is None):
            return None
        value = str(value).strip()
        if(value.lstrip('-').isdigit()):
            index = int(value)
            if(not 0 <= index < len(self.layers) - 1):
                self._fail(f"{name} = {index}: the layers below the sticky air have IDs 0 to {len(self.layers) - 2}.")
            return index
        indices = layer_indices(self.layers, value)
        if(not indices):
            self._fail(f"{name} = {value}: no unit labelled '{value}' in {self.svg_file}.")
        return max(indices)

    def _insert_basalt(self, index):
        '''Insert a basalt layer of zero thickness at index of the stack (below the layer at index); return its index.'''
        basalt = MandyocLayer('basalt', BASALT, Nx=self.Nx, title='Basalt', color='xkcd:dark grey blue',
                              rho=2900.0, C=1.0, H=9.0e-11)
        basalt.interface = np.copy(self.layers[index].interface) #base of basalt = base of the layer above it
        basalt.base_depth = basalt.interface.max()
        self.layers.insert(index, basalt)
        return index

    def _read_settings(self, settings):
        where = f"Setting of {self.svg_file}"
        def setting(name, default, quantity=None, unit=None, scale=1.0):
            '''
            Setting (argument > SVG > default). A text with a unit of tapioca.UNITS[quantity] (without unit: unit) is
            converted to the base unit of the quantity and multiplied by scale (e.g. yr -> Myr).
            '''
            value = self._setting(name, settings.get(name), default)
            if(quantity is not None and isinstance(value, str)):
                value = parse_quantity(value, quantity, where, name, default_unit=unit) * scale
            return value

        self.regime = self._check_choice('regime', setting('regime', 'rifting'), REGIMES)
        self.velocity = self._check_number('velocity', setting('velocity', 1.0, 'velocity', 'cm/yr')) #cm/yr
        self.moving_boundary = self._check_choice('moving_boundary', setting('moving_boundary', 'right'),
                                                  ['left', 'right', 'symmetric'])
        time_max = setting('time_max', None, 'time', 'Myr', 1.0e-6) #Myr
        if(time_max is not None):
            time_max = self._check_number('time_max', time_max, 0, strict=True)
        self.quiescence_time = self._check_number('quiescence_time', setting('quiescence_time', 0.0, 'time', 'Myr', 1.0e-6), 0) #Myr

        self.magmatism = self._check_on_off('magmatism', setting('magmatism', 'off'))
        self.magmatism_extraction = self._check_on_off('magmatism_extraction', setting('magmatism_extraction', 'off'))
        self.magmatism_Ts0 = self._check_number('magmatism_Ts0', setting('magmatism_Ts0', 1080.0, 'temperature', 'oC')) #oC
        self._check_experiment(self.regime, self.moving_boundary, self.source, self.magmatism, self.magmatism_extraction)

        preset = setting('preset_initial_temperature', None)
        if(preset is not None): #a relative path in the SVG is relative to the folder of the SVG
            base = os.path.dirname(os.path.abspath(self.svg_file)) if(self.source['preset_initial_temperature'] == 'SVG') else None
            preset = self._check_path('preset_initial_temperature', preset, base)
        self.preset_initial_temperature = preset

        self.temperature_solver = self._check_choice('temperature_solver', setting('temperature_solver', '2D'),
                                                     ['2D', 'columns'], aliases={'1D': 'columns'})

        #Velocity and time:
        #   rifting: extension, the moving boundary moves outwards with velocity
        #   compression: convergence, the moving boundary pushes the lithosphere inwards with velocity
        #   accordion (as in M300): the boundary velocity of the rifting regime is scaled in time by scale_bcv.txt
        #       0 - ti_convergence: extension
        #       ti_convergence - ti_quiescence: convergence (velocity x -1)
        #       ti_quiescence - tf_quiescence: quiescence (velocity x 0.01)
        #       tf_quiescence - time_max: extension again (dt_rifting2)
        #   quiescence: no velocity imposed at the boundaries (velocity_from_ascii = False)
        self.velocity_from_ascii = (self.regime != 'quiescence')
        self.variable_bcv = (self.regime == 'accordion')
        self.dt_quiescence = self.quiescence_time
        if(self.variable_bcv):
            self.tf_quiescence = self.ti_quiescence + self.dt_quiescence
            if(time_max is not None): #time_max sets the duration of the second extension
                self.dt_rifting2 = time_max - self.tf_quiescence
                if(self.dt_rifting2 <= 0):
                    self._fail(f"time_max = {time_max:g} Myr ends before the second extension starts ({self.tf_quiescence:g} Myr).")
            self.time_max = (self.tf_quiescence + self.dt_rifting2) * 1.0e6
        elif(time_max is not None):
            self.time_max = time_max * 1.0e6
        elif(self.regime == 'quiescence'):
            self.time_max = self.dt_quiescence * 1.0e6 if(self.dt_quiescence > 0) else 1000.0e6
        else: #rifting and compression
            self.time_max = 60.0e6

        self._print_verbose(f"Numerical experiment: {self.regime} ({self.source['regime']}), velocity = {self.velocity:g} cm/yr "
                    f"({self.source['velocity']}), moving boundary: {self.moving_boundary} ({self.source['moving_boundary']})")
        self._print_verbose(f"Magmatism: {self.magmatism} ({self.source['magmatism']}), extraction: {self.magmatism_extraction} "
                    f"({self.source['magmatism_extraction']}), Ts0 = {self.magmatism_Ts0:g} oC ({self.source['magmatism_Ts0']})")
        if(self.preset_initial_temperature is not None):
            self._print_verbose(f"Initial temperature from the preliminary scenario {self.preset_initial_temperature} "
                        f"({self.source['preset_initial_temperature']})")

    ###########################################################################
    # Profiles for the plots
    ###########################################################################

    ## I think its a very complex function for a simple task, all data is stored inside the fields dataset
    ## All these profiles could be obtained from the: MandyocBuilder.DTree.fields.ds
    def _find_profiles(self):
        '''
        Temperature and YSE profiles of the plots: one for each lithospheric domain (e.g. oceanic,
        non-cratonic, cratonic), i.e. each region with the same lithosphere thickness and the same unit at its
        base over at least min_profile_width (outside the weak seed). Each profile is at the centre of the widest
        region and is labelled with the title of the lithospheric mantle unit above the asthenosphere.
        Return a list of (column, label, color).
        '''
        layers, Nx = self.layers, self.Nx
        litho_km = np.round(self.thickness_litho / 1.0e3, 1)
        has_seed = np.zeros(Nx, dtype=bool)
        for i in layer_indices(layers, self.seed_label):
            has_seed |= slot_thickness(layers, i, slice(None)) > 0
        base_unit = np.array([next((layers[i].layer_label for i in range(self.idx_ast + 1, len(layers))
                                    if slot_thickness(layers, i, ix) > 0), '') for ix in range(Nx)])
        litho_key = np.array([f"{litho_km[ix]}|{base_unit[ix]}" if(not has_seed[ix]) else '' for ix in range(Nx)])
        edges = np.r_[0, np.nonzero(litho_key[1:] != litho_key[:-1])[0] + 1, Nx]
        regions = {} #(lithosphere thickness, unit at its base): (first column, last column + 1)
        for begin, end in zip(edges[:-1], edges[1:]):
            key = (litho_km[begin], base_unit[begin])
            if(litho_key[begin] != '' and (end - begin) * self.dx >= self.min_profile_width * self.Lx):
                if((key not in regions) or (end - begin > np.diff(regions[key])[0])):
                    regions[key] = (begin, end)

        profiles = []
        colors = self.profile_colors[-len(regions):] if(len(regions) <= len(self.profile_colors)) else self.profile_colors
        for k, ((thickness, unit), (begin, end)) in enumerate(sorted(regions.items())):
            ix = (begin + end - 1) // 2
            i_mantle = next((i for i in range(self.idx_ast + 1, len(layers)) if slot_thickness(layers, i, ix) > 0), self.idx_ast + 1)
            name = layers[i_mantle].title.replace(' lith. mantle', '').replace(' lithospheric mantle', '')
            profiles.append((ix, f"{name} ({thickness:.0f} km)", colors[k % len(colors)]))
        if(len(profiles) == 0):
            profiles = [(self.idx_profile, f"x = {self.x_profile/1.0e3:.0f} km", 'r')]
        return profiles

    ###########################################################################
    # Initial temperature
    ###########################################################################

    def build_temperature(self):
        '''
        Initial temperature (Nz, Nx) built with a TemperatureFieldBuilder (self.temperature_builder), whose
        lithosphere goes from the surface (base of the sticky air) to the top of the asthenosphere (LAB):

        - without preset_initial_temperature: steady-state conduction in the lithosphere, with the
          conductivity and the radiogenic heat of each layer, 0 oC at the surface and TP at the LAB, and the
          mantle adiabat below the LAB (TemperatureFieldBuilder.apply_linear_temperature). temperature_solver
          '2D' (default) solves the whole lithosphere at once (solve_heat_2D), so lateral conduction follows the
          geometry of the layers; 'columns' solves each column alone (solve_heat_1D);
        - with preset_initial_temperature: the final thermal state of the preliminary scenario adapted to the
          lithosphere of each column.

        The field is also the temperature of the builder (self.builder.DTree.fields.temperature).
        '''
        self.temperature_builder = TemperatureFieldBuilder(self.builder, LAB_interface=self.idx_ast,
                                                           surface_interface=self.surface_id, c_cap=self.ccapacity)
        if(self.preset_initial_temperature is None):
            self.temperature_builder.apply_linear_temperature(t_pot=self.TP)
            if(self.temperature_solver == '2D'):
                self.temperature_builder.solve_heat_2D()
            else:
                self.temperature_builder.solve_heat_1D()
        else:
            try:
                depth, T_preset, sa_preset, litho_preset, self.preset_source = read_preset_temperature(self.preset_initial_temperature)
            except FileNotFoundError as err:
                raise SVGInputError(f"preset_initial_temperature: {err}")
            #Without the interfaces.txt of the preliminary scenario, its lithosphere is the most common one of this scenario
            if(litho_preset is None):
                values, counts = np.unique(np.round(self.thickness_litho / 1.0e3), return_counts=True)
                litho_preset = values[np.argmax(counts)] * 1.0e3
                print(f"Warning: interfaces.txt of the preliminary scenario not found; its lithosphere is assumed to be "
                      f"{litho_preset/1.0e3:.0f} km thick (the most common in {self.svg_file}).")
            self.thickness_litho_preset = litho_preset
            self.temperature_builder.t_pot = self.TP
            #sticky air and lithosphere of each column from the surface and the LAB of the temperature builder
            tb = self.temperature_builder
            tb.temperature.values = adapt_preset_temperature(depth, T_preset, sa_preset, litho_preset, -tb.z,
                                                             -tb.surface_interface.values, tb.lithosphere_thickness.values,
                                                             self.thickness_mantle)
            self._print_verbose(f"Preset of initial temperature: {self.preset_source}")
            self._print_verbose(f"\tsticky air of the preliminary scenario: {sa_preset/1.0e3:.0f} km; reference lithosphere: "
                        f"{litho_preset/1.0e3:.0f} km (other lithosphere thicknesses are adapted)")
        self.builder.create_temperature_field(self.temperature_builder)
        self.T = self.builder.DTree.fields.temperature.values #(Nz, Nx)
        return self.T

    ###########################################################################
    # Velocity boundary conditions
    ###########################################################################

    def build_velocity(self):
        '''
        Initial velocity (Nz, Nx) with the boundary conditions of the regime (None in the quiescence regime),
        set by VelocityFieldBuilder.set_kinematic_regime (self.velocity_builder): constant velocity down to h_v_const
        below the surface (reference thickness + 20 km, where the lithosphere is thickest in rifting and accordion,
        as in M300, and at the moving boundary in compression) and the return flow below it, so that the integral of
        the inflow is zero. Rifting and accordion are extensions; in compression the return flow leaves through the
        moving boundary and the opposite boundary is fixed. With moving_boundary = symmetric (rifting and accordion),
        both boundaries move outwards with velocity/2, each with its own return flow. In accordion, scale_bcv.txt
        changes the velocity in time.
        The velocity fields are those of the builder (self.builder.DTree.fields.vx and vz).
        '''
        if(not self.velocity_from_ascii):
            self.VX, self.VZ = None, None
            return None

        idx_bc = 0 if(self.moving_boundary == 'left') else self.Nx - 1 #column of the moving boundary
        thickness_ref_x = self.thickness_litho if(self.h_v_reference == 'lithosphere') else self.thickness_mantle
        thickness_ref_bc = thickness_ref_x[idx_bc] if(self.regime == 'compression') else thickness_ref_x.max()
        self.h_v_const = thickness_ref_bc + 20.0e3 #thickness with constant velocity

        v = VelocityFieldBuilder(self.builder, units='cm/y')
        v.set_kinematic_regime('compression' if(self.regime == 'compression') else 'extension', self.velocity,
                              moving_boundary=self.moving_boundary, h_const=self.h_v_const, air_gap=self.fac_air,
                              surface=self.surface_id)
        self.builder.create_velocity_field(v)
        self.velocity_builder = v

        self.VX = self.builder.DTree.fields.vx.values
        self.VZ = self.builder.DTree.fields.vz.values
        return self.VX, self.VZ

    def scale_bcv(self):
        '''Content of scale_bcv.txt in the accordion regime (Mandyoc multiplies the boundary velocity by each factor at the given time).'''
        if(not self.variable_bcv):
            return None
        if(self.dt_quiescence > 0):
            changes = [(self.ti_convergence, -1.0), (self.ti_quiescence, 0.01), (self.tf_quiescence, -100.0),
                       (self.time_max/1.0e6, 0.01)]
        else: #no quiescence: convergence followed directly by extension
            changes = [(self.ti_convergence, -1.0), (self.ti_quiescence, -1.0), (self.time_max/1.0e6, 0.01)]
        return f"{len(changes)}\n" + "".join(f"{t} {factor}\n" for t, factor in changes)

    def parameters(self, **params):
        '''
        Parameters of param.txt (dict): DEFAULT_PARAMS (tapioca.inputs.svg._variables) and those of the scenario
        (time_max, velocity and magmatism settings). The grid (nx, nz, lx, lz) and n_interfaces are set by
        MandyocBuilder.set_params. Any parameter of PARAMS_TYPES can be changed with a keyword argument,
        e.g. parameters(particles_per_element=100, rheology_model=9).
        '''
        unknown = [key for key in params if key not in PARAMS_TYPES]
        if(unknown):
            raise SVGInputError(f"Unknown parameter(s) of param.txt: {', '.join(sorted(unknown))}")
        return {**DEFAULT_PARAMS,
                'time_max': self.time_max,
                'variable_bcv': self.variable_bcv,
                'velocity_from_ascii': self.velocity_from_ascii,
                'magmatism': self.magmatism,
                'magmatism_extraction': self.magmatism_extraction,
                'magmatism_Ts0': self.magmatism_Ts0,
                'magmatic_layer': self.magmatic_layer,
                **params}

    ###########################################################################
    # Writing the files
    ###########################################################################

    def write_scripts(self, output_dir='.', system=None, ncores=None, mandyoc_options='',
                      petsc_dir=None, petsc_arch=None, mandyoc=None):
        '''
        Write the scripts of the scenario in output_dir:
          - run_linux.sh or run_mac.sh: runs Mandyoc with mpirun of PETSc and then moves the outputs to one
            folder per variable with 0_organize_outputs.sh (on the Mac, caffeinate keeps it awake);
          - 0_organize_outputs.sh: copy of src/postrunning/0_organize_outputs.sh;
          - zipper.sh: zips the scenario in <folder>.zip, skipping the files already in the zip.

        Parameters
        ----------
        system : str, optional
            'linux' or 'mac'. Default: the system that runs tapioca.
        ncores : int, optional
            Number of MPI processes. Default: number of cores of this computer.
        mandyoc_options : str, optional
            Command line options of Mandyoc (e.g. '-seed 0,5,8 -strain_seed 0.0,1.0,1.0').
        petsc_dir, petsc_arch, mandyoc : str, optional
            PETSc folder and architecture (mpirun is $PETSC_DIR/$PETSC_ARCH/bin/mpirun) and the Mandyoc
            executable. Default: PETSC_DIR and PETSC_ARCH of the environment, or $HOME/opt/petsc and the
            only arch-* folder in it with mpirun (arch-label-optimized if there is none or more than one),
            and $HOME/opt/mandyoc/bin/mandyoc. Values set in the environment when the script runs override them.
        '''
        this_system = {'Linux': 'linux', 'Darwin': 'mac'}.get(platform.system())
        if(system is None):
            if(this_system is None):
                raise SVGInputError(f"Run scripts are written for {' or '.join(SYSTEMS)}; give system "
                                    f"(this computer is {platform.system()}).")
            system = this_system
        if(system not in SYSTEMS):
            raise SVGInputError(f"Unknown system '{system}'; use {' or '.join(SYSTEMS)}.")
        if(not os.path.isfile(ORGANIZE_SCRIPT)):
            raise SVGInputError(f"0_organize_outputs.sh not found in {os.path.dirname(os.path.abspath(ORGANIZE_SCRIPT))}.")

        if(petsc_dir is None):
            petsc_dir = os.environ.get('PETSC_DIR', '$HOME/opt/petsc')
        if(petsc_arch is None):
            petsc_arch = os.environ.get('PETSC_ARCH')
        if(petsc_arch is None):
            #the arch of $HOME/opt/petsc, only if the script is for this computer
            archs = glob.glob(os.path.join(os.path.expanduser('~'), 'opt', 'petsc', 'arch-*', 'bin', 'mpirun'))
            petsc_arch = (archs[0].split(os.sep)[-3] if(system == this_system and len(archs) == 1)
                          else 'arch-label-optimized')
        if(mandyoc is None):
            mandyoc = '$HOME/opt/mandyoc/bin/mandyoc'
        if(ncores is None):
            ncores = os.cpu_count() or 1

        os.makedirs(output_dir, exist_ok=True)
        path = lambda fname: os.path.join(output_dir, fname)

        run = RUN_SCRIPT
        for key, value in [('PETSC_DIR', petsc_dir), ('PETSC_ARCH', petsc_arch), ('MANDYOC', mandyoc),
                           ('NCORES', int(ncores)), ('MANDYOC_OPTIONS', mandyoc_options),
                           ('LAUNCHER', RUN_LAUNCHER[system])]:
            run = run.replace(f'@{key}@', str(value))
        scripts = {f'run_{system}.sh': run, 'zipper.sh': ZIP_SCRIPT}
        with open(ORGANIZE_SCRIPT) as f:
            scripts['0_organize_outputs.sh'] = f.read()

        for fname, text in scripts.items():
            with open(path(fname), 'w') as f:
                f.write(text)
            os.chmod(path(fname), 0o755)
        if(self.verbose):
            print(f"Scripts: {', '.join(scripts)} ({system}, {int(ncores)} cores)")

    def write_all(self, output_dir='.', plots=True, copy_svg=True, scripts=True, system=None, ncores=None,
                  mandyoc_options='', **params):
        '''
        Write the input files of Mandyoc in output_dir: interfaces.txt, param.txt, input_temperature_0.txt,
        input_velocity_0.txt, vel_bc.txt, velz_bc.txt and scale_bcv.txt (accordion), the scenario infos
        (infos_<folder>.txt), a copy of the SVG and the figures numerical_setup.png and initial_temperature.png.
        With scripts, also the run script of Mandyoc for system ('linux' or 'mac'), 0_organize_outputs.sh and
        zipper.sh (see write_scripts).
        Keyword arguments change parameters of param.txt (see parameters).
        '''
        os.makedirs(output_dir, exist_ok=True)
        path = lambda fname: os.path.join(output_dir, fname)

        if(copy_svg):
            target = path(os.path.basename(self.svg_file))
            if(not os.path.exists(target) or not os.path.samefile(self.svg_file, target)):
                shutil.copy(self.svg_file, target)

        builder = self.builder
        builder.path = output_dir
        builder.export_interfaces(export='none')
        builder.set_params(**self.parameters(**params))
        builder.export_params()

        if(self.T is None):
            self.build_temperature()
        builder.export_field('temperature', header="T1\nT2\nT3\nT4")

        if(self.velocity_from_ascii):
            if(self.VX is None):
                self.build_velocity()
            #fields of the builder, which may have been changed by script after build_velocity
            self.VX, self.VZ = builder.DTree.fields.vx.values, builder.DTree.fields.vz.values
            X = self.X
            v0 = self.VX[(X == 0)]
            vf = self.VX[(X == self.Lx)]
            np.savetxt(path("vel_bc.txt"), np.array([v0, vf]).T)
            np.savetxt(path("velz_bc.txt"), self.VZ[(self.z == 0)].T)
            builder.export_field('velocity', header="v1\nv2\nv3\nv4")
            if(self.variable_bcv):
                with open(path("scale_bcv.txt"), "w") as f:
                    f.write(self.scale_bcv())

        if(scripts):
            self.write_scripts(output_dir, system=system, ncores=ncores, mandyoc_options=mandyoc_options)

        if(plots):
            self.plot_temperature(path('initial_temperature.png'))
            self.plot_setup(path('numerical_setup.png'))

        name = os.path.basename(os.path.abspath(output_dir))
        infos = self.infos(name)
        np.savetxt(path(f'infos_{name}.txt'), infos, fmt="%s")
        if(self.verbose):
            for line in infos:
                print(line)

    ###########################################################################
    # Scenario infos
    ###########################################################################

    ## This is a useful function that could be inside the MandyocBuilder class
    def infos(self, name=''):
        '''Description of the scenario (list of lines).'''
        layers, ix = self.layers, self.idx_profile
        lines = ['SCENARIO INFOS:', ' ', f'Name: {name}',
                 f"Built with tapioca {__version__} from {os.path.basename(self.svg_file)}: SVG {self.svg_version or '(not declared)'}, "
                 f"saved by Inkscape {self.inkscape_version.split()[0] if self.inkscape_version else '(not declared)'} "
                 f"(tested: SVG {SVG_VERSION}, Inkscape {INKSCAPE_VERSION})",
                 ' ', 'Domain parameters:',
                 f"\tLx: {self.Lx*1.0e-3} km", f"\tLz: {self.Lz*1.0e-3} km", f"\tNx: {self.Nx}", f"\tNz: {self.Nz}",
                 f"Resolution dx x dz: {1.0e-3*self.dx} x {1.0e-3*self.dz} km2", ' ',
                 f'Time limit: {self.time_max/1.0e6} Myr']
        if(self.variable_bcv):
            lines += [f"Extension: 0 - {self.ti_convergence} Myr",
                      f"Convergence: {self.ti_convergence} - {self.ti_quiescence} Myr",
                      f"Time of quiescence: {self.dt_quiescence} Myr ({self.ti_quiescence} - {self.tf_quiescence} Myr)",
                      f"Second extension: {self.dt_rifting2} Myr ({self.tf_quiescence} - {self.time_max/1.0e6} Myr)", ' ']
        lines.append(f'Layers (bottom to top) from {self.svg_file}, thickness at x = {self.x_profile/1.0e3:.0f} km:')
        for i, layer in enumerate(layers):
            lines.append(f"\t{i}: {layer.layer_label} ({layer.rheology_name}): h = {slot_thickness(layers, i, ix)/1.0e3:.1f} km, "
                         f"C = {layer.params['C']}, rho = {layer.params['rho']} kg/m3, H = {layer.params['H']:.4g} W/kg")
        lines.append(f"\tlithosphere: {self.thickness_litho[ix]/1.0e3:.1f} km")
        lines.append(f"\tmagmatic layer: {self.magmatic_layer} ({layers[self.magmatic_layer].layer_label})")
        lines.append("Temperature and YSE profiles in numerical_setup.png:")
        for ixp, label, color in self.profiles:
            lines.append(f"\t{label}: x = {self.x[ixp]/1.0e3:.0f} km")
        lines += [' ', f"Numerical experiment: {self.regime}", f"Velocity field: {self.velocity_from_ascii}"]
        if(self.velocity_from_ascii and self.VX is not None):
            lines.append(f"initial velocity: {self.velocity:g} cm/yr")
            lines.append(f"\tmoving boundary: {self.moving_boundary}, constant velocity down to {self.h_v_const/1.0e3:.0f} km below the surface")
        lines += [f"Variable velocity field: {self.variable_bcv}", 'Initial temperature field setup:',
                  f"\tPreset of initial temperature field: {self.preset_initial_temperature is not None}"]
        if(self.preset_initial_temperature is None):
            lines += [f"\tSteady-state conduction in the lithosphere ({self.temperature_solver} solver)",
                      f"\tMantle potential temperature (at the LAB): {self.TP} oC"]
        elif(hasattr(self, 'preset_source')):
            lines += [f"\tPreliminary scenario: {self.preset_source}",
                      f"\tReference lithosphere of the preliminary scenario: {self.thickness_litho_preset/1.0e3:.0f} km"]
        lines += [' ', f'magmatism: {self.magmatism}, extraction: {self.magmatism_extraction}, Ts0: {self.magmatism_Ts0:g} oC']
        lines.append('Settings (argument, SVG or default): ' + ', '.join(f"{k}: {v}" for k, v in self.source.items()))
        return lines

    ###########################################################################
    # Plots
    ###########################################################################

    def yield_strength_envelope(self, ix, sr=1.0e-15):
        '''YSE [Pa] along the vertical profile at column ix (strain rate sr [1/s]), from MandyocBuilder.calculate_YSE.'''
        if(self.T is None):
            self.build_temperature()
        fields = self.builder.DTree.fields
        if('YSE' not in fields or fields.YSE.attrs.get('strain_rate') != sr):
            self.builder.calculate_YSE(sr=sr)
        return self.builder.DTree.fields.YSE.isel(x=ix).values

    ## It's a good function to be inside the TemperatureFieldBuilder
    ## This function should get **kwargs to set the plot
    ## TemperatureFieldBuilder has all interfaces encapsulated
    def plot_temperature(self, fname='initial_temperature.png'):
        '''Initial temperature field with the 600, 1000 and 1300 oC isotherms and the layer interfaces.'''
        T = self.T if(self.T is not None) else self.build_temperature()
        X, Z, x = self.X, self.Z, self.x
        sa = self.thickness_sa[self.idx_profile]

        fig, ax0 = plt.subplots(nrows=1, ncols=1, constrained_layout=True, figsize=(16, 6))
        im = ax0.contourf(X / 1.0e3, (-Z + sa) / 1.0e3, T, levels=np.arange(0, np.max(T) + 100, 100), cmap='inferno')
        isotherms = ax0.contour(X / 1.0e3, (-Z + sa) / 1.0e3, T, levels=[600, 1000, 1300], colors='w', linewidths=1.0, linestyles='--')
        ax0.clabel(isotherms, fmt='%d °C', fontsize=10)
        for layer in self.layers[1:]:
            ax0.plot(x / 1.0e3, (-layer.interface + sa) / 1.0e3, color='k', lw=0.6)
        ax0.set_xlim(0, self.Lx / 1.0e3)
        ax0.set_ylim(-self.Lz/1000 + sa/1000, 0 + sa/1000)
        ax0.tick_params(labelsize=12)
        ax0.set_xlabel("Distance [km]", fontsize=14)
        ax0.set_ylabel("Depth [km]", fontsize=14)
        ax0.set_title(f"Initial temperature field ({os.path.basename(self.svg_file)})", fontsize=14)
        cbar = fig.colorbar(im, orientation='vertical', ax=ax0, pad=0.01)
        cbar.set_label("Temperature [°C]", fontsize=12)
        cbar.ax.tick_params(labelsize=10)
        fig.savefig(fname, dpi=300)
        plt.close(fig)

    def plot_setup(self, fname='numerical_setup.png'):
        '''
        Numerical setup: layers with the colours of the SVG, legend with the properties, weak seed,
        temperature and YSE profiles of each lithospheric domain, friction angle and velocity boundary conditions.
        '''
        T = self.T if(self.T is not None) else self.build_temperature()
        if(self.velocity_from_ascii and self.VX is None):
            self.build_velocity()
        with plt.rc_context({'xtick.labelsize': 18, 'ytick.labelsize': 18}):
            self._plot_setup(fname, T)

    def _plot_setup(self, fname, T):
        layers, x, z, Lx, Lz, Nx = self.layers, self.x, self.z, self.Lx, self.Lz, self.Nx
        sa = self.thickness_sa[self.idx_profile]
        dx, dz = self.dx, self.dz

        fig, axs = plt.subplots(1, 1, figsize=(14, 6))
        ylimplot = [-Lz/1000 + sa/1000, 0 + sa/1000]

        for i, layer in enumerate(layers):
            top = layers[i+1].interface if i < len(layers)-1 else np.zeros(Nx)
            axs.fill_between(x/1000, (-layer.interface + sa)/1000, (-top + sa)/1000, color=layer.color, lw=0)

        axs.set_xlim(0, Lx/1000)
        axs.set_xticks([])
        axs.set_yticks([])
        axs.set_ylim(ylimplot)
        axs.set_xlabel(f"Distance = {Lx/1000:.0f} km\nResolution = {dx/1000:.0f} km", fontsize=16)
        axs.set_ylabel(f"Depth = {Lz/1000:.0f} km\nResolution = {dz/1000:.0f} km", fontsize=16)

        #legend: one entry per label, from top to bottom; thickness h:
        #   bodies (units that do not reach both sides of the domain, e.g. a basin): maximum thickness (depocentre)
        #   bottom unit (asthenosphere): range of thickness along x
        #   other layers: most common thickness along x (their reference thickness)
        legend_layers = {}
        for layer in layers[::-1]:
            if(layer.layer_label != self.seed_label):
                legend_layers.setdefault(layer.layer_label, layer)

        def legend_thickness(label):
            th = layer_thickness(layers, label, slice(None)) / 1.0E3
            if(th[0] <= 0 or th[-1] <= 0):
                return f"{th.max():.0f} km"
            if(label == layers[0].layer_label and np.round(th.min()) != np.round(th.max())):
                return f"{th.min():.0f}-{th.max():.0f} km"
            values, counts = np.unique(np.round(th[th > 0], 1), return_counts=True)
            return f"{values[np.argmax(counts)]:.0f} km" if(len(values)) else "0 km"

        legend_elements = []
        for label, layer in legend_layers.items():
            legend_label = (f"{layer.title}\n{layer.params['C']:g} x {layer.rheology_name.replace('_', ' ')}\n"
                            + fr"$\rho$ = {layer.params['rho']:.0f} kg/m³" + f"\n$h$ = {legend_thickness(label)}")
            legend_elements.append(Patch(facecolor=layer.color, edgecolor='black', label=legend_label))
        fig.subplots_adjust(bottom=0.25)
        axs.legend(handles=legend_elements, ncol=len(legend_elements), loc='lower center', bbox_to_anchor=(0.5, -0.45),
                   frameon=False, title='Layers properties', title_fontsize=15, fontsize=8, columnspacing=1.0)


        ## Idea: get the weak seed from the MandyocBuilder.DTree.interfaces, where weakening_seed > 0
        #weak seed (one annotation per body, i.e. per group of contiguous columns)
        idx_seed = layer_indices(layers, self.seed_label)
        if(len(idx_seed) > 0):
            seed_bases = np.array([layers[i].interface for i in idx_seed])
            seed_tops = np.array([layers[i+1].interface if i < len(layers)-1 else np.zeros(Nx) for i in idx_seed])
            seed_present = seed_bases > seed_tops
            seed_top = np.where(seed_present, seed_tops, np.inf).min(axis=0)
            seed_base = np.where(seed_present, seed_bases, -np.inf).max(axis=0)
            seed_cols = seed_present.any(axis=0)
            seed_title = layers[idx_seed[0]].title.lower()
            edges = np.diff(np.r_[0, seed_cols.astype(int), 0])
            for begin, end in zip(np.nonzero(edges == 1)[0], np.nonzero(edges == -1)[0]):
                cols = slice(begin, end)
                x_seed = np.mean(x[cols])
                width_seed = (end - begin) * dx
                thickness_seed = np.max(seed_base[cols] - seed_top[cols])
                axs.text(x_seed/1000, -(np.max(seed_base[cols]) - sa)/1000 - 2,
                         f'{seed_title}\n{width_seed/1000:.0f}x{thickness_seed/1000:.0f} km²',
                         color='k', fontsize=12, ha='center', va='top')

        #temperature profiles
        axt = axs.inset_axes((0.205, 0, 0.18, 1))
        for ix, label, color in self.profiles:
            axt.plot(T[:, ix], (-z + sa) / 1.0e3, "-", color=color, label=label)
            axs.plot(x[ix] / 1.0e3, sa / 1.0e3 * 0.15, 'v', color=color, ms=6) #position of the profile
        axt.legend(loc='lower left', fontsize=8, framealpha=0.9)
        axt.grid(visible=True, axis='x', which='both', ls='--', color='red', alpha=0.3)
        axt.set_ylim(ylimplot)
        axt.set_yticks([])
        axt.set_xticks(np.linspace(0, 1800, 7))
        axt.patch.set_alpha(0)
        axt.xaxis.set_ticks_position('top')
        axt.tick_params(labelsize=8)
        axt.tick_params(axis='x', colors='k')
        axt.spines['left'].set_visible(False)
        axt.spines['right'].set_visible(False)
        axt.set_title('Temperature [°C]', color='k')

        #yield strength envelopes
        axsg = axs.inset_axes((0.605, 0, 0.13, 1))
        for ix, label, color in self.profiles:
            axsg.plot(self.yield_strength_envelope(ix)/1e9, -(z - sa)/1e3, '-', color=color, label=label)
        axsg.grid(visible=True, axis='x', which='both', ls='--', color='gray', alpha=0.8)
        axsg.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        axsg.set_xlim(-0.1, 1.1)
        axsg.set_ylim(ylimplot)
        axsg.set_yticks([])
        axsg.patch.set_alpha(0)
        axsg.xaxis.set_ticks_position('top')
        axsg.spines['left'].set_visible(False)
        axsg.spines['right'].set_visible(False)
        axsg.set_title(r'$\sigma_{YSE}$ [GPa]')
        axsg.tick_params(labelsize=8)

        #effective friction angle
        axsf = axs.inset_axes((0.800, 0, 0.10, 0.30))
        axsf.patch.set_alpha(0.4)
        fisize = 10
        axsf.plot(np.array([0, 0.25, 0.75, 1.0]), np.array([15, 15, 2, 2]), 'k-')
        axsf.set_xlim([0.10, 0.9])
        axsf.set_ylim([0, 17])
        axsf.set_xticks([0.25, 0.5, 0.75])
        axsf.set_xticklabels([0.05, ' ', 1.05])
        axsf.tick_params('x', top=True, labeltop=False)
        axsf.xaxis.set_ticks_position('bottom')
        axsf.set_xlabel(r"$\varepsilon$", fontsize=fisize)
        axsf.set_yticks([2, 15])
        axsf.set_yticklabels(['2°', '15°'])
        axsf.set_ylabel(r'$\Phi_{\mathrm{eff}}$', fontsize=fisize)
        axsf.tick_params(labelsize=fisize)
        axsfC = axsf.twinx()
        axsfC.patch.set_visible(False) #the background of the twin axes would hide the curve
        axsfC.set_ylim([0, 17])
        axsfC.set_yticks([2, 15])
        axsfC.set_yticklabels([4, 20])
        axsfC.set_ylabel('Cohesion [MPa]', fontsize=fisize, rotation=90)
        axsfC.tick_params(labelsize=fisize)

        #velocity boundary conditions: insets centred on the edges, so v = 0 is on the boundary
        if(self.velocity_from_ascii):
            scale = 100.0 * SEC_PER_YEAR #m/s to cm/yr
            fac = 0.96
            for column, origin, title in ((-1, fac, f'v = {self.velocity} [cm/y]'), (0, -(1-fac), None)):
                vmax = np.round(max(abs(self.VX[:, column] * scale)), 0)
                vmax = vmax if(vmax > 0) else 1.0
                axv = axs.inset_axes((origin, 0, (1-fac)*2, 1))
                #only the hatch is drawn: the hatch takes the edge color (with alpha=0 it is also invisible in matplotlib >= 3.11)
                axv.fill_betweenx((-z + sa) / 1.0e3, self.VX[:, column] * scale, 0, facecolor='none', edgecolor='k', linewidth=0, hatch='---')
                axv.set_ylim(ylimplot)
                axv.set_yticks([])
                axv.set_xticks([])
                axv.patch.set_alpha(0)
                if(title):
                    axv.set_title(title, fontsize=8)
                axv.set_xlim([-vmax, vmax])
                for side in ['left', 'right', 'bottom'] + (['top'] if(title is None) else []):
                    axv.spines[side].set_visible(False)
            #the hatches of the left boundary go beyond the domain where VX < 0 (clipped by the inset): the y label
            #(its right side is at the given x) is 0.5% of the width to the left of them
            vmin_left = max(-self.VX[:, 0].min() * scale, 0.0)
            hatch_left = (1 - fac) * min(vmin_left / vmax, 1.0)
            axs.yaxis.set_label_coords(-hatch_left - 0.005, 0.5)

        fig.savefig(fname, bbox_inches="tight", dpi=300)
        plt.close(fig)
