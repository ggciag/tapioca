import numpy as np

_variables = ["density",
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
"surface",
"Phi",
"dPhi",
"X_depletion"]

_varsTypes = { # Outputs
              'density': np.float64,
              'pressure': np.float64,
              'heat': np.float64,
              'thermal_difussivity': np.float64,
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


_seg_per_year = 60*60*24*365.25 # seconds per year

cm = 1/2.54 # cm per inch -> convert plotting to centimeters by *cm 