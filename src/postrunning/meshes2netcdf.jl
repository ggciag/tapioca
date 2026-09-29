# Write by Joao Bueno - Jun. 2026
# This script converts text files from Mandyoc outputs into netcdf4 files.
# Outputs must be organized in folders (check `0_organize_outputs.sh`)


using NCDatasets
using Glob
using CSV
using Printf
using DataFrames
using Base.Threads
using StatsBase

global AIR_DENSITY_THRESHOLD = -999
global LITHOLOGY_DATATYPE = Int8
global VARIABLES = ["density", "viscosity", "pressure", "strain","strain_rate","temperature","velocity","heat"]
global CHUNKS = 5
global dfllevel = 7 # compression level 1-9

global UNITS =Dict{String,String}(
    "x"=>"m",
    "y"=>"m",
    "z"=>"m",
    "time"=>"Myr",

    "density"=>"kg/m3",
    "viscosity"=>"Pa.s",
    "pressure"=>"Pa",
    "strain"=>"dimensionless",
    "strain_rate"=>"1/s",
    "temperature"=>"deg C",
    "velocity"=>"m/s",
    "surface"=>"m",
    "heat"=>"W/kg",
    "thermal_diffusivity"=>"m2/s",
    "Phi"=>"dimensionless",
    "dPhi"=>"1/s",  # Need to confirm
    "X_depletion"=>"dimensionless",
    "lithology"=>"ID",
)

global DTYPES =Dict{String,DataType}(
    "x"=>Float32,
    "y"=>Float32,
    "z"=>Float32,
    "time"=>Float32,

    "density"=>Float32,
    "viscosity"=>Float64,
    "pressure"=>Float64,
    "strain"=>Float64,
    "strain_rate"=>Float64,
    "temperature"=>Float32,
    "velocity"=>Float64,
    "surface"=>Float32,
    "heat"=>Float64,
    "thermal_diffusivity"=>Float64,
    "Phi"=>Float64,
    "dPhi"=>Float64,  # Need to confirm
    "X_depletion"=>Float64,
    "lithology"=>LITHOLOGY_DATATYPE,
)

# Useful for non dimensional scenarios
global SCALE_FACTOR = Dict{String,Number}(


)

struct mesh2D 
    Nx::Int
    Nz::Int
    Lx::Float64
    Lz::Float64
end

struct mesh3D
    Nx::Int
    Ny::Int
    Nz::Int
    Lx::Float64
    Ly::Float64
    Lz::Float64
end

struct MandyocScenario
    dims::Int
    steps::Vector{Int32}
    times::Vector{Float64}
    thick_air::Float32
    units::Dict{String,String}
    datatypes::Dict{String,DataType}
end

function read_param(fpath::String="param.txt")::Dict{String,String}
    param_dict = Dict{String,String}()
    open(fpath, "r") do file
        for line in eachline(file)
            line=strip(line)
            if isempty(line) || startswith(line, "#")
                continue
            end
            line = split(line, "#")[1]
            line = replace(line, " " => "")
            key_value = split(line, "=")
            if length(key_value) == 2
                param_dict[lowercase(key_value[1])] = key_value[2]
            end
        end
    end
    return param_dict
end

function read_data(var::String, step::Integer, mesh::mesh2D, vartype::Type; veloc::Bool=false, surface::Bool=false)
    
    file = "$(var)/$(var)_$(step).txt"
    Nx,Nz = mesh.Nx, mesh.Nz
    
    skipto::Int = 0
    if surface skipto=0 else skipto=3 end
    C = CSV.File(file, header=false, comment="P", skipto=skipto,types=vartype)|>CSV.Tables.matrix
    
    
    if veloc 
        C[abs.(C) .< 1e-50] .= 0
        vx = transpose(reshape(C[1:2:end], (Nx, Nz)))
        vz = transpose(reshape(C[2:2:end], (Nx, Nz)))
        R = (vx,vz)
        return R
    
    elseif surface
        return (C[1:end,1],C[1:end,2]) #sx, sy
    
    else
        C[abs.(C) .< 1e-50] .= 0
        R = transpose(reshape(C[1:end], (Nx, Nz)))
        return R
    end
end

function get_all_steps()::Vector{Int32}
    pattern = joinpath("time", "time_*.txt")
    files = glob(pattern)
    
    steps = Int32[parse(Int32, split(splitext(basename(f))[1], '_')[end]) for f in files]
    return sort(steps)
end

function read_time(step::Integer)::Float64
    time::Float64=0.0
    open(joinpath("time", "time_$step.txt")) do file
        line = readline(file)
        line = split(line,"   ")[2]
        time = parse(Float64,line)/1e6 #Myr
    return time
    end
end

# Create the netcdf for the original mesh
function create_nc(variable::String,scen::MandyocScenario, mesh::mesh2D)

    Nx, Nz = mesh.Nx, mesh.Nz
    Lx, Lz = mesh.Lx, mesh.Lz
    times, steps = scen.times, scen.steps

    nc_fname = "$(variable).nc"
    
    if variable == "lithology"
        Nx = (Nx-1)*5 + 1
        Nz = (Nz-1)*5 + 1
    end

    x_coords = DTYPES["x"].(range(0.0f0, Lx, length=Nx))
    z_coords = DTYPES["z"].(range(-Lz, 0.0f0, length=Nz))
    
    num_steps = length(steps)
    vardtype = DTYPES[variable]

    Dataset(nc_fname,"c") do ds #criar o arquivo nc

        defDim(ds,"time",num_steps)
        defVar(ds,"time", times, ("time",),
        attrib=Dict("units"=>UNITS["time"],"long_name"=>"time","axis"=>"T"),
                                                                deflatelevel=dfllevel, shuffle=true)
        
        defVar(ds,"step", Int32, ("time",), attrib=Dict("units"=>"","long_name"=>"step"),
                deflatelevel=dfllevel, shuffle=true)
        ds["step"][:] = steps

        if variable == "surface"
            sx_sample,_ = read_data("surface", 0, mesh, DTYPES["surface"], veloc=false, surface=true)
            Nxs = size(sx_sample)[1]
            surf_x = DTYPES["x"].( range(0.0f0, Lx, length=Nxs) )
            defDim(ds,"x",Nxs)
            defVar(ds,"x",surf_x,("x",),attrib=Dict("units"=>UNITS["x"],"long_name"=>"x","axis" => "X"), 
                                                                deflatelevel=dfllevel, shuffle=true)
            
            defVar(ds, variable, vardtype, ("x", "time"),
                attrib=Dict("long_name"=>variable, "units"=>get(UNITS, variable, "-")),
                deflatelevel=dfllevel, shuffle=true)
        else
            defDim(ds,"x",Nx)
            defDim(ds,"z",Nz)

            defVar(ds,"x",x_coords,("x",),attrib=Dict("units"=>UNITS["x"],"long_name"=>"x","axis" => "X"),
                                                                    deflatelevel=dfllevel, shuffle=true,
                                                                    )
            defVar(ds,"z",z_coords,("z",),attrib=Dict("units"=>UNITS["z"],"long_name"=>"z","axis"=>"Z"),
                                                                    deflatelevel=dfllevel, shuffle=true)
            
            if variable == "velocity"
                defVar(ds,"vx",vardtype,("x","z","time"),attrib=Dict("units"=>UNITS[variable],
                                                                    "long_name"=>"vx"),
                                                                    deflatelevel=dfllevel, shuffle=true)
                
                defVar(ds,"vz",vardtype,("x","z","time"),attrib=Dict("units"=>UNITS[variable],
                                                                    "long_name"=>"vz",),
                                                                    deflatelevel=dfllevel, shuffle=true)
            
            else 
                defVar(ds, variable, vardtype, ("x", "z", "time"),
                attrib=Dict(
                    "long_name"=>variable,
                    "units"=>get(UNITS, variable, "-")
                ),
                deflatelevel=dfllevel, shuffle=true)
            end
        end
    end
end

function converter(variable::String, scen::MandyocScenario, mesh::mesh2D)
    nc_fname = "$variable.nc"
    Nx = mesh.Nx
    Nz = mesh.Nz
    Lx = mesh.Lx
    Lz = mesh.Lz
    steps = scen.steps
    times = scen.times
    vardtype = DTYPES[variable]
    num_steps = length(steps)

    veloc = (variable == "velocity")
    surface = (variable == "surface")
    
    # Create the NC file
    create_nc(variable,scen,mesh)
    
    # Iterate chunks -> open NC file
    chunk_size = ceil(Int, num_steps / CHUNKS) # Define amount of steps per chunk
    step_chunks = collect(Iterators.partition(1:num_steps, chunk_size)) # Partitions the indices safely

    for (chunk_idx, indices) in enumerate(step_chunks)
        sub_steps = steps[indices]
        len_chunk = length(sub_steps) 
        println("Running chunk $chunk_idx/$CHUNKS - Steps $(sub_steps[1]) to $(sub_steps[end])")
        
        # Allocate buffers just for this chunk
        local buffer_vx, buffer_vz, buffer_var, buffer_surf
        if surface
            sx_sample, _ = read_data("surface", sub_steps[1], mesh, DTYPES["surface"], veloc=false, surface=true)
            buffer_surf = zeros(DTYPES["surface"], size(sx_sample)[1], len_chunk)
        elseif veloc
            buffer_vx = zeros(DTYPES["velocity"], Nx, Nz, len_chunk)
            buffer_vz = zeros(DTYPES["velocity"], Nx, Nz, len_chunk)
        else
            buffer_var = zeros(vardtype, Nx, Nz, len_chunk)
        end

        progress_counter = Threads.Atomic{Int}(0)
        start_time = time()
        
        @threads for i in eachindex(sub_steps)
            step = sub_steps[i]
            data = read_data(variable, step, mesh, vardtype, veloc=veloc, surface=surface)
            
            if data !== nothing
                if veloc
                    dens = read_data("density", step, mesh, DTYPES["density"], veloc=false, surface=false)
                    vx, vz = data
                    vx[dens .< AIR_DENSITY_THRESHOLD] .= 0
                    vz[dens .< AIR_DENSITY_THRESHOLD] .= 0
                    buffer_vx[:, :, i] = vx'
                    buffer_vz[:, :, i] = vz'
                elseif surface
                    _, sy = data
                    buffer_surf[:, i] = sy
                else
                    dens = read_data("density", step, mesh, DTYPES["density"], veloc=false, surface=false)
                    data[dens .< AIR_DENSITY_THRESHOLD] .= 0
                    buffer_var[:, :, i] = data'
                end
            else
                @warn "No data found for $(variable)_$(step).txt at step $step"
            end
            
            Threads.atomic_add!(progress_counter, 1)
            if progress_counter[] % 10 == 0
                speed = (time() - start_time) / progress_counter[]
                active_ram = round(Base.gc_live_bytes() / 1000^2, digits=1)
                ram_peak = round(Sys.maxrss() / 1000^2, digits=1)
                @info "[$variable] Progress: $(progress_counter[])/$len_chunk | Speed: $(round(speed, digits=2))s/step | Active RAM: $(active_ram)MB | RAM peak: $(ram_peak)MB"
            end
        end
        
        # Append chunk data to NC file
        Dataset(nc_fname, "a") do ds
            if veloc
                ds["vx"][:, :, indices] = buffer_vx
                ds["vz"][:, :, indices] = buffer_vz
            elseif surface
                ds["surface"][:, indices] = buffer_surf
            else
                ds[variable][:, :, indices] = buffer_var
            end
        end
        
        # Free chunk memory
        buffer_vx = buffer_vz = buffer_var = buffer_surf = nothing
        GC.gc()
    end
    println("Saved to $nc_fname\n---------------")
end

# ======= Functions to convert Lithology =======

function read_litho_file(fpath::String)
    # Read lithology file and return x, z, lith
    conf = CSV.File(fpath, 
                    header=false, 
                    comment="P", 
                    delim=' ', 
                    types=[Int32, Int32, Int8], # x, z, lith_id
                    silencewarnings=true)
                    
    return Vector(conf.Column1), Vector(conf.Column2), Vector(conf.Column3)
end

function replace_negatives_with_neighbors!(mat::Matrix)
    result = mat
    rows::Int16, cols::Int16 = size(mat)
    negative_indices = findall(result .< 0)

    @inbounds for idx in negative_indices
        i::Int16, j::Int16 = idx[1], idx[2]
        found = false
        @inbounds for (di, dj) in [(-1, 0), (1, 0), (0, -1), (0, 1)]
            ni, nj = i + di, j + dj
            if (1 <= ni <= rows && 1 <= nj <= cols) && mat[ni, nj] >= 0
                result[i, j] = mat[ni, nj]
                found = true
                break
            end
        end
    end
    return result
end

function convert_litho_to_nc(scen::MandyocScenario, mesh::mesh2D, cores::Integer)
    nc_fname = "lithology.nc"
    Nx, Nz = mesh.Nx, mesh.Nz
    Lx, Lz = mesh.Lx, mesh.Lz
    steps = scen.steps
    times = scen.times
    ncores = cores
    num_steps = length(steps)
    
    # Upscaled mesh for lithology
    Nxl = (Nx - 1) * 5 + 1
    Nzl = (Nz - 1) * 5 + 1
    
    create_nc("lithology", scen, mesh)

    chunk_size = ceil(Int, num_steps / CHUNKS)
    step_chunks = collect(Iterators.partition(1:num_steps, chunk_size))

    start_time_global = time()

    for (chunk_idx, indices) in enumerate(step_chunks)
        sub_steps = steps[indices]
        len_chunk = length(sub_steps) 
        println("Running chunk $chunk_idx/$CHUNKS - Steps $(sub_steps[1]) to $(sub_steps[end])")

        buffer = zeros(LITHOLOGY_DATATYPE, Nxl, Nzl, len_chunk)
        total_steps = length(steps)

        progress_counter = Threads.Atomic{Int}(0)
        chunk_start_time = time()

        @threads for i in eachindex(sub_steps)
            step = sub_steps[i]
             # LITHOLOGY_DATATYPE can't be unssigned if -1 is the fill value!
            litho_mesh = fill(LITHOLOGY_DATATYPE(-1), Nzl, Nxl)
            
            for core in 0:(ncores-1)
                fpath = joinpath("lithos", "litho_$(step)_$core.txt")
                if isfile(fpath)
                    x_core, z_core, litho_core = read_litho_file(fpath)
                    for k in eachindex(x_core)
                        litho_mesh[z_core[k] + 1, x_core[k] + 1] = litho_core[k]
                    end
                end
            end

            replace_negatives_with_neighbors!(litho_mesh)

            buffer[:, :, i] = reverse(litho_mesh, dims=1)'
            
            Threads.atomic_add!(progress_counter, 1)
            if progress_counter[] % 10 == 0
                speed = (time() - chunk_start_time) / progress_counter[]
                ram_peak = round(Sys.maxrss() / 1000^2, digits=1)
                active_ram = round(Base.gc_live_bytes() / 1000^2, digits=1)
                @info "[lithology] Progress: $(progress_counter[])/$len_chunk | Speed: $(round(speed, digits=2))s/step | Active RAM: $(active_ram)MB | RAM peak: $(ram_peak)MB"
            end
        end 

        Dataset(nc_fname, "a") do ds
            ds["lithology"][:, :, indices] = buffer
        end
        
        # cleaning the ram
        buffer = nothing
        GC.gc()
    end
    
    total_elapsed = time() - start_time_global
    @info "Finished! Total time: $(round(total_elapsed / 60, digits=2)) minutes."
end

function build_scenario(params::Dict)

    dims = parse(Int, get(params, "dimensions", "2")) # if the model is 2d or 3d

    # Change data types in case of non dimensional scenarios 
    if get(params,"iterative","direct") == "iterative" || get(params,"nondimensionalization","False") == "True" || dims == 3
        for v in keys(DTYPES) DTYPES[v] = Float64 end
        global AIR_DENSITY_THRESHOLD = -1
        println("You did run a non dimensional model, all datatypes changed to Float64.")
    end

    Nx = parse(Int, params["nx"]) # elements in x
    Nz = parse(Int, params["nz"]) # elements in z
    Lx = parse(DTYPES["x"], params["lx"]) # x length
    Lz = parse(DTYPES["z"], params["lz"]) # z length
    thick_air::DTYPES["z"] = 40.0f3 # m

    if dims == 3 
        Ny = parse(Int, params["ny"])
        Ly = parse(Int, params["ly"])
        mesh = mesh3D(Nx,Ny,Nz,Lx,Ly,Lz)
    else
        mesh = mesh2D(Nx,Nz,Lx,Lz)
    end

    # Finding all steps
    steps = get_all_steps()
    times = DTYPES["time"][read_time(step) for step in steps]
    CSV.write("times.csv", DataFrame(step=steps, time_myr=times))
    println("$(length(steps)) time steps were found, from $(times[1]) Myr [$(steps[1])] up to $(times[end]) Myr [$(steps[end])].")

    return MandyocScenario(dims, steps, times, thick_air, UNITS, DTYPES), mesh
end

function main()
    
data_dir = ARGS[end] # Scenario directory
cd(data_dir)

# Basic parameters
params = read_param("param.txt")

if (get(params,"sp_surface_tracking", "False") == "True") || (get(params,"sp_surface_processes", "False") == "True")
    push!(VARIABLES, "surface")
    println("Surface was tracked.")
end

if get(params, "magmatism", "off") == "on"
    push!(VARIABLES, "Phi")
    push!(VARIABLES, "dPhi")
    push!(VARIABLES, "X_depletion")
    println("Magmatism was on.")
end

if get(params, "export_thermal_diffusivity", "False") == "True"
    push!(VARIABLES, "thermal_diffusivity")
    println("Thermal diffusivity (kappa) was exported.")
end

scen, mesh = build_scenario(params)

for var in unique(VARIABLES)
    println("Converting: $var")
    converter(var,scen,mesh)
end

println("All variables were converted to NetCDF4")

if get(params, "export_lithology", "False") == "True"
    println("Lithology grid was exported.")
    ncores::Int = size(glob(joinpath("lithos","litho_0_*.txt")))[1]
    println("$ncores cores were used in this model.")
    convert_litho_to_nc(scen,mesh, ncores)
    println("Lithology was converted to NetCDF4")
end


println("Finished")
println("Variables converted: $(join(VARIABLES,"; "))")
println("Compression level: $(dfllevel)")
println("-"^20)
println("Types:")
println(join(DTYPES,";\n"))
println("-"^20)

end

main()
