# Post running scripts

These scripts contain code to facilitate the data formatting that TapIOca expects. Most of these scripts were made to handle with .txt files. Note that .txt files will be become deprecated in future versions of Mandyoc.

The common routine is (1) to organise text file outputs and (2) to convert them into netcdf files.

An example:

```
bash 0_organize_outputs.sh
julia -t THREADS meshes2netcdf.jl SCENARIO_PATH
```
in which `THREADS` is the number of desired threads, and `SCENARIO_PATH` is the path to the modelled scenario folder.