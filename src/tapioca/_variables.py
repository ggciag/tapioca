import numpy as np

_variables_list = ["density",
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


seg_per_year = 60*60*24*365.25 # seconds per year

cm = 1/2.54 # cm per inch -> convert plotting to centimeters by *cm 


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