'''
Materials of Mandyoc scenarios: layers (MandyocLayer) and their rheologies (LithologicalUnit).

The rheologies are instances of LithologicalUnit built from tapioca.DEFAULT_MATERIALS
(WET_OLIVINE, DRY_OLIVINE, WET_QUARTZ, ...), so adding a material to that dictionary makes it
available here and in the SVG descriptions.
'''

import numpy as np

from ..._variables import DEFAULT_MATERIAL, DEFAULT_MATERIALS

#CLASSES

## I think the MandyocLayer could be only a interface (xr.DataArray) in the MandyocBuilder.DTree.interfaces (xr.Dataset)
class MandyocLayer:
    def __init__(self, layer_label, rheology, Nx, base_depth=0.0, title=None, color="#808080",
                 reference_thickness=None, **params):
        """
        A layer of the stack of an SVG scenario: its material (params), its interface (depth of its base along
        x) and the attributes of the drawing (title, color, reference thickness h).

        Parameters
        ----------
        layer_label: str
            Name of the layer
        rheology: LithologicalUnit or str
            Rheology of the layer, e.g. WET_OLIVINE or its name ("wet olivine", see get_rheology)
        Nx: int
            Number of points in x direction
        base_depth: float
            Depth of the layer base [m] (the interface is flat until it is set)
        title: str, optional
            Name of the layer in the legend of the figures (Title of the shape in the SVG). Default: layer_label.
        color: str, optional
            Color of the layer in the figures (fill of the shape in the SVG). Default: grey.
        reference_thickness: float, optional
            Reference thickness h of the unit [m], given in the SVG description (None if not given).
        params: keyword arguments
            Parameters of the layer with the keys of tapioca.DEFAULT_MATERIAL (rho, C, H, A, n, Q, V, k,
            weakening_seed, cohesion_min, ...); they replace the ones of the rheology.

        Attributes
        ----------
        params: dict
            The 13 parameters of interfaces.txt (those of the rheology, replaced by the ones given).
        given: dict
            The parameters given to the layer (e.g. in the description of the SVG).
        """
        self.layer_label = layer_label
        self.rheology = rheology if(isinstance(rheology, LithologicalUnit)) else get_rheology(rheology)
        self.rheology_name = self.rheology.name
        unknown = set(params) - set(DEFAULT_MATERIAL)
        if unknown:
            raise ValueError(f"Unknown parameter(s) of {layer_label}: {', '.join(sorted(unknown))}. Use {', '.join(DEFAULT_MATERIAL)}.")
        self.given = dict(params)
        self.params = {**self.rheology.params, **params}
        self.base_depth = base_depth
        self.Nx = Nx
        self.title = title if title is not None else layer_label
        self.color = color
        self.reference_thickness = reference_thickness
        self.interface = np.ones(Nx)*base_depth

    def material(self):
        '''Parameters of the layer in interfaces.txt (the keys of tapioca.DEFAULT_MATERIAL), as given to MandyocBuilder.create_interface.'''
        return dict(self.params)

class LithologicalUnit:
    """
    Generic material: a name and the parameters of interfaces.txt (the keys of tapioca.DEFAULT_MATERIAL:
    rho, C, H, A, n, Q, V, k, cohesion_min, ...). Parameters not given take the DEFAULT_MATERIAL values.

    Example: LithologicalUnit("wet_olivine", A=1.393e-14, n=3.0, Q=429.0e3, V=14.5e-6)

    name: str
        Name of the material, e.g. wet_olivine
    params: keyword arguments
        Parameters of the material
    """
    def __init__(self, name, **params):
        unknown = set(params) - set(DEFAULT_MATERIAL)
        if unknown:
            raise ValueError(f"Unknown parameter(s) of {name}: {', '.join(sorted(unknown))}. Use {', '.join(DEFAULT_MATERIAL)}.")
        self.name = name
        self.params = {**DEFAULT_MATERIAL, **params}

    def __repr__(self):
        return f"LithologicalUnit({self.name!r}, {', '.join(f'{k}={v}' for k, v in self.params.items())})"


#Default rheologies, one instance of LithologicalUnit per material of DEFAULT_MATERIALS (WET_OLIVINE, AIR, ...)
RHEOLOGIES = {key.lower(): LithologicalUnit(key.lower(), **params) for key, params in DEFAULT_MATERIALS.items()}
globals().update({unit.name.upper(): unit for unit in RHEOLOGIES.values()})

__all__ = ["MandyocLayer", "LithologicalUnit", "RHEOLOGIES", "get_rheology"] + [name.upper() for name in RHEOLOGIES]

def get_rheology(name):
    '''
    Return the rheology (LithologicalUnit) from a free-form name.
    "wet olivine", "wet_olivine", "Wet Olivine" and "WetOlivine" all return WET_OLIVINE.
    '''
    key = name.replace(" ", "").replace("_", "").replace("-", "").lower()
    for rheology in RHEOLOGIES.values():
        if rheology.name.replace("_", "") == key:
            return rheology
    options = ", ".join(RHEOLOGIES)
    raise ValueError(f"Unknown rheology '{name}'. Options: {options}")
