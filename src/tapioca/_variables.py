import numpy as np

__all__ = ["VARIABLES_LIST","VARS_TYPES",
           "SEC_PER_YEAR","CM",
           "INTERFACES_PARAMETERS","DEFAULT_MATERIAL","MATERIAL_PARAMETERS","PARAMETERS_UNITS"
           "DEFAULT_MATERIALS"]

VARIABLES_LIST:list = ["density",
"heat",
"strain",
"lithology",
"strain_rate",
"surface",
"pressure",
"temperature",
"thermal_diffusivity",
"velocity",
"viscosity",
"Phi",
"dPhi",
"X_depletion"]
"""list of str: The standard list of expected Mandyoc output variables."""

VARS_TYPES:dict = { # Outputs
              'density': np.float64,
              'pressure': np.float64,
              'heat': np.float64,
              'thermal_diffusivity': np.float64,
              'surface': np.float64,
              'lithology': np.int8, # int 8 is enough for 256 lithologies, and it saves memory
              'viscosity': np.float64, 
              'velocity': np.float64,
              'strain': np.float64,
              'strain_rate':np.float64,

              # Post processing
              'deviatoric_stress':np.float64, # considering \\tau in Pa
              
              #Dimensions
              'time': np.float64,
              'x': np.float64,
              'z': np.float64,
              'id': np.int64 # for particles
}
"""dict: Dictionary with datatypes (value) for each variable (key)."""

SEC_PER_YEAR:float = 60*60*24*365.25 # seconds per year
"""float: Conversion factor of seconds in a year."""

CM:float = 1/2.54 # cm per inch -> convert plotting to centimeters by *cm 
"""float: Conversion factor of inches to centimeters. Useful to change matplotlib measures."""

INTERFACES_PARAMETERS:dict = {
    "C": "scale_factor",
    "rho": "density",
    "H": "radiogenic_heat",
    "A": "pre-exponential_scale_factor",
    "n": "power_law_exponent",
    "Q": "activation_energy",
    "V": "activation_volume",
    "k": "thermal_diffusivity",
    "weakening_seed": "weak_seed_strain",
    "cohesion_min": "cohesion_min",
    "cohesion_max": "cohesion_max",
    "friction_angle_min": "friction_angle_min",
    "friction_angle_max": "friction_angle_max"
}
"""dict: Dictionary containing the interfaces.txt variables and their corresponding aliases."""

DEFAULT_MATERIAL:dict = {
    "rho": None,  
    "C": 1.0,
    "H": 0.0,
    "A": 1e-18,
    "n": 1,
    "Q": 0.0,
    "V": 0.0,
    "k": 1.0e-6,
    "weakening_seed": -1,
    "cohesion_min": 4e6,
    "cohesion_max": 20e6,
    "friction_angle_min": 2,
    "friction_angle_max": 15
}
"""dict: Dictionary containing default properties for the interface material (linear rheology)"""

MATERIAL_PARAMETERS:list = list(DEFAULT_MATERIAL.keys())

PARAMETERS_UNITS:dict = {
    "rho" : 'kg/m3',
    "C" : 'dimensionless',
    "H" : 'W/kg',
    "A" : 'Pa^(-n)/s',
    "n" : 'dimensionless',
    "Q" : 'J/mol',
    "V" : 'm3/mol',
    "k" : 'm2/s',
    "weakening_seed" : 'dimensionless',
    "cohesion_min" : 'Pa',
    "cohesion_max" : 'Pa',
    "friction_angle_min": 'degrees',
    "friction_angle_max": 'degrees',
}
"""dict: Dictionary containing the units for each interface parameters."""

DEFAULT_MATERIALS:dict = {
    "WET_OLIVINE": { # Karato and Wu (1993)
        "C": 1.0,
        "rho": 3378.0,
        "H": 0.0,
        "A": 1.393e-14,
        "n": 3.0,
        "Q": 429000.0,
        "V": 1.450e-05,
        "k": 1.0e-6,
        "weakening_seed": -1.0,
        "cohesion_min": 4000000.0,
        "cohesion_max": 20000000.0,
        "friction_angle_min": 2.0,
        "friction_angle_max": 15.0
    },
    "DRY_OLIVINE": { # Karato and Wu (1993)
        "C": 1.0,
        "rho": 3354.0,
        "H": 9.000e-12,
        "A": 2.417e-15,
        "n": 3.5,
        "Q": 540000.0,
        "V": 2.450e-05,
        "k": 1.0e-6,
        "weakening_seed": -1.0,
        "cohesion_min": 4000000.0,
        "cohesion_max": 20000000.0,
        "friction_angle_min": 2.0,
        "friction_angle_max": 15.0
    },
    "WET_QUARTZ": {  # Glean and Tullis (1995)
        "C": 1.0,
        "rho": 2700.0,
        "H": 4.630e-10,
        "A": 8.574e-28,
        "n": 4.0,
        "Q": 222000.0,
        "V": 0.0,
        "k": 1.0e-6,
        "weakening_seed": -1.0,
        "cohesion_min": 4000000.0,
        "cohesion_max": 20000000.0,
        "friction_angle_min": 5.0,
        "friction_angle_max": 15.0
    },
    "WET_ANORTHITE": { # Rybacki and Dresen (2000); Andrés-Martínez et al. (2019)
        "C": 1.0,
        "rho": 2850.0,
        "H": 7.123e-11,        # 0.041 uW/m^3
        "A": 8.913e-22,     # 10^-21.05
        "n": 3.0,
        "Q": 356000.0,      # 356 kJ/mol
        "V": 0.0,
        "k": 1.15e-6,
        "weakening_seed": -1.0,
        "cohesion_min": 4000000.0,
        "cohesion_max": 25000000.0,
        "friction_angle_min": 2.0,
        "friction_angle_max": 25.0
    },
    "LINEAR_SALT": { # Massimi et al. (2007); Pichel et al. (2022)
        "C": 1.0,
        "rho": 2200.0,
        "H": 0.0,
        "A": 2.0e-19,
        "n": 4.0,
        "Q": 0.0,
        "V": 0.0,
        "k": 3.5e-6,
        "weakening_seed": -1.0,
        "cohesion_min": 4000000.0,
        "cohesion_max": 20000000.0,
        "friction_angle_min": 2.0,
        "friction_angle_max": 15.0
    },
    "AIR": { 
        "C": 0.1,
        "rho": 1.0,
        "H": 0.0,
        "A": 1.0e-18,
        "n": 1.0,
        "Q": 0.0,
        "V": 0.0,
        "k": 1.0e-5,
        "weakening_seed": -1.0,
        "cohesion_min": 100000.0,
        "cohesion_max": 1000000.0,
        "friction_angle_min": 2.0,
        "friction_angle_max": 5.0
    }
}
'''dict: Dictionary with the parameters of some "classic" rheologies. The rheology is the key and the values is the parameters.'''


# Useful data types descriptions:

'''
int8
range: -128 to 127
size: 1 byte
precision: none

int16
range: -32,768 to 32,767
size: 2 bytes
precision: none

int32
range: -2,147,483,648 to 2,147,483,647 (-2^31 to 2^31 - 1)
size: 4 bytes
precision: none

int64
range: -9.22e18 to 9.22e18 (-2^63 to 2^63 - 1)
size: 8 bytes
precision: none

uint8
range: 0 to 255
size: 1 byte
precision: none

uint16
range: 0 to 65,535
size: 2 bytes
precision: none

float16
range: -65,504 to +65,504
size: 2 bytes
precision: ~3 decimal digits

float32
range: -3.4e38 to +3.4e38
size: 4 bytes
precision: ~7 decimal digits

float64
range: -1.79e308 to +1.79e308
size: 8 bytes
precision: ~15 to 17 decimal digits
'''