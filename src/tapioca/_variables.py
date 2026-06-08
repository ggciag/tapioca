import numpy as np
#Vars
_varsTypes = { # Outputs
              'density': np.float16, #smaller float / float
              'pressure': np.float64,
              'heat': np.float64,
              'thermal_difussivity': np.float64,
              'surface': np.float32,
              'lithology': np.int8, #
              'viscosity': np.float64, #
              'velocity': np.float64,
              'strain': np.float64,
              'strain_rate':np.float64,


              # Post processing
              'deviatoric_stress':np.float32, # considers \\tau in MPa
              
              #Dimensions
              'time': np.float32,
              'x': np.float32,
              'z': np.float32,
              'id': np.int64 # for particles
}