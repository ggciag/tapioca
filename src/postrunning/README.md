# Post running scripts

These scripts contain code to facilitate the data formatting that TapIOca expects. Most of these scripts were made to handle with .txt files. Note that .txt files will be become deprecated in future versions of Mandyoc.

The common routine is (1) to organise text file outputs and (2) to convert them into netcdf files.

An example:

```
bash 0_organize_outputs.sh
julia -t THREADS meshes2netcdf.jl SCENARIO_PATH
```
in which `THREADS` is the number of desired threads, and `SCENARIO_PATH` is the path to the modelled scenario folder. 
To create the dataset files, you can filter variables in the sticky air layer by setting the variable `AIR_DENSITY_THRESHOLD` with the air density used (usually 1-10 kg/m³). This will set all densities lower than this threshold to zero.
The user can specify a number of chuncks to process the outputs using `CHUNKS`. It is useful for fine grids that can consume dozens of RAM memory. Additionaly, the user is free to set the data types for each variable using the dict `DTYPES`. 

The amount of necessary RAM can be estimated using the velocity and lithology datasets. Velocity is the most memory-intensive standard variable because it requires two 8 bytes arrays (vx and vz), and reading it also requires loading the density array (4 bytes) into memory simultaneously. 

- Buffer Memory: $\text{Nx} \times \text{Nz} \times \text{Steps\_per\_Chunk} \times 16 \text{ bytes}$
- Thread Overhead: $\text{Active\_Threads} \times \text{Nx} \times \text{Nz} \times 20 \text{ bytes}$ (Loading vx, vz, and density into the thread)
- Total: $\text{Max RAM (Velocity)} = \frac{\text{Nx} \times \text{Nz}}{1024^2} \times (16 \times \text{Steps\_per\_Chunk} + 20 \times \text{Active\_Threads})$ MB.

Lithology uses 1 byte by default (Int8) which is very light, but the output is upscaled by a factor of 5 in both axis. This means the lithology grid has roughly 25 times more elements than your original mesh. 
- Upscaled Elements: $\text{Nxl} \times \text{Nzl} \approx 25 \times (\text{Nx} \times \text{Nz})$
- Buffer Memory: $\text{Nxl} \times \text{Nzl} \times \text{Steps\_per\_Chunk} \times 1 \text{ byte}$
- Thread Overhead: $\text{Active\_Threads} \times \text{Nxl} \times \text{Nzl} \times 1 \text{ byte}$ (Plus a small CSV parsing overhead per thread)
- Total: $\text{Max RAM (Lithology)} = \frac{\text{Nxl} \times \text{Nzl}}{1024^2} \times (\text{Steps\_per\_Chunk} + \text{Active\_Threads})$ MB.

In which $\text{Steps\_per\_Chunk} = \text{Total\_Steps} / \text{CHUNKS}$. A pratical amount of necessary RAM for each chunk is: $1.25 \times \max{R_{\text{litho}},R_{\text{veloc}}}$.