# Post-Processing

Post-processing objects and classes are in root directory. 

## Mandyoc Scenario Class

Creating a `MandyocScen` object. The source code is in `scenClasses.py`.

```{eval-rst}
.. currentmodule:: tapioca.scenClasses

.. autosummary::
   :toctree: ../generated/

   MandyocScen

```

### Methods

Methods of the `MandyocScen` class.

```{eval-rst}

.. currentmodule:: tapioca.scenClasses

.. autosummary::
   :toctree: ../generated/

   MandyocScen.correctZcoord

```

#### Internal methods

```{eval-rst}

.. currentmodule:: tapioca.scenClasses

.. autosummary::
   :toctree: ../generated/

   MandyocScen._load_spatial_var
   MandyocScen._load_particles


```

## Post-Processing Accessors and Functions

Post-processing acessors and functions are in `post_processing.py`. 

The `scenPostProcessing` class is a DataTree accessor for calculating common post-processing tasks and routines. 
This class registers under the ``postproc`` namespace for ``MandyocScen`` class providing methods to compute physical properties directly on the DataTree structure.

```{eval-rst}
.. currentmodule:: tapioca.post_processing

.. autosummary::
   :toctree: ../generated/

   scenPostProcessing.DeviatoricStressTensor
```

---

Example:
```python
scen = MS('example-scen', name='Example',
                   variables=['density','temperature','viscosity','velocity'])
scen.DTree.postproc.DeviatoricStressTensor()
```