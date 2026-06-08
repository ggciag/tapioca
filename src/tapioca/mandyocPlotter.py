
import json
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from .scenClasses import *
from ._variables import *

class mandyocPlotter: 
    '''
    [Under construction] A class to plot mandyoc data.
    '''
    
    def __init__(self, mandyocScen, 
                 cmap_dir='/media/jobueno/STOV/scripts/salt_cmap.json'):
        
        return self
    
        
    def read_cmap(self, cmap_dir):
        json_dir = cmap_dir
        file = open(json_dir,"r")
        cmap_json = json.load(file)
        self.cmap_metadata = cmap_json['metadata']
        self.cmap = LinearSegmentedColormap.from_list(name=self.cmap_metadata['name'], 
                                                      colors=cmap_json['colors'], 
                                                      N=cmap_json['metadata']['N_layers'])
        #plt.register_cmap(name=cmap_metadata['name'], cmap=cmap_json['colors'])
        file.close()
        return True
    
    class snapshot(): #sub classe de Plotter para plotar os campos/cenarios de 1 step selecionado
        
        def __init__(self, tselect):
        
            
            return self