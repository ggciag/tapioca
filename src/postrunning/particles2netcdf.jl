# Written by Joao Bueno - Dec. 2025
# Updated by Joao Bueno - Sep. 2026

# This script converts steps text files from Mandyoc into a netcdf4 file following the CF Convention.
# Outputs must be organized in sub-folders (check `0_organize_outputs.sh`)

# ======= CONSTANTS =======

const LITHOLOGY_DATATYPE = Int8         # set it to Int16 if you have more than 127 lithologies
const CHUNKS = 2                        # set it to optimize ram usage and/or processing time
const DEFLATE_LEVEL = 5                 # set it to optimize the file compression (from 1 to 9)
const SELECT_LITHOLOGY = Int8.(0:127)   # set it to only select particles with the given lithologies

# it possible to select particles within a domain using the variables "SELECT_X" and "SELECT_Z" in the main()

global DTYPES =Dict{String,DataType}(
    "x"=>Float32,
    "y"=>Float32,
    "z"=>Float32,
    "id"=>Int32,
    "lithology"=>LITHOLOGY_DATATYPE
)

# ======= FUNCTIONS AND STRUCTS =======

using NCDatasets
using Glob
using CSV
using Printf
using DataFrames
using Base.Threads

function read_param(fpath::String="param.txt")::Dict{String,String}
    param_dict = Dict{String,String}()
    open(fpath, "r") do file
        for line in eachline(file)
            line = strip(line)
            if isempty(line) || startswith(line, "#") continue end
            key_value = split(replace(split(line, "#")[1], " " => ""), "=")
            if length(key_value) == 2
                param_dict[lowercase(key_value[1])] = key_value[2]
            end
        end
    end
    return param_dict
end

function get_all_steps()::Vector{Int32}
    files = glob(joinpath("time", "time_*.txt"))
    steps = Int32[parse(Int32, split(splitext(basename(f))[1], '_')[end]) for f in files]
    return sort(steps)
end

function read_time(step::Integer)::Float64
    time = 0.0
    open(joinpath("time", "time_$step.txt")) do file
        line = split(readline(file), "   ")[2]
        time = parse(Float64, line) / 1e6 # Myr
    end
    return time
end

function get_ncores()::Int
    return length(glob(joinpath("steps", "step_0_*.txt"))) 
end

# ======= getting the reference IDs =======

function extract_reference_ids(ref_step::Int, ncores::Int)
    println("Loading the reference time step: ($ref_step)")
    dfs = Vector{DataFrame}(undef, ncores)
    
    @threads for core in 0:(ncores-1)
        fpath = joinpath("steps", "step_$(ref_step)_$(core).txt")
        if isfile(fpath)
            dfs[core+1] = CSV.read(fpath, DataFrame, 
                                   header=false, delim=' ', ignorerepeated=true,
                                   select=[1, 2, 3, 4], 
                                   types=[DTYPES["x"], DTYPES["z"], DTYPES["id"], DTYPES["lithology"], Float64], 
                                   silencewarnings=true)
        end
    end
    
    df_ref = vcat(dfs...)
    rename!(df_ref, [:x, :z, :id, :lithology])
    
    # Filtering
    filter!(row -> 
        row.x >= SELECT_X[1] && row.x <= SELECT_X[2] &&
        row.z >= SELECT_Z[1] && row.z <= SELECT_Z[2] &&
        row.lithology in SELECT_LITHOLOGY, 
    df_ref)

    sort!(df_ref, :id)
    target_ids = df_ref.id
    
    # Map array
    max_id = maximum(target_ids)
    id_map = zeros(DTYPES["id"], max_id + 1)
    for (i, id) in enumerate(target_ids)
        id_map[id + 1] = i
    end

    println("Total of tracked particles: $(length(target_ids))")
    return target_ids, id_map
end

function create_particles_nc(nc_fname::String, target_ids::Vector{Int32}, steps::Vector{Int32}, times::Vector{Float64}, ref_step::Int)
    num_particles = length(target_ids)
    num_steps = length(steps)

    Dataset(nc_fname, "c") do ds
        defDim(ds, "id", num_particles)
        defDim(ds, "time", num_steps)

        # 1D Coordinates
        defVar(ds, "id", target_ids, ("ID",), attrib=Dict("long_name"=>"particle_id",
                                                          "cf_role"   => "trajectory_id"), 
                                        deflatelevel=DEFLATE_LEVEL) 

        defVar(ds, "time", times, ("time",), attrib=Dict("units"=>"myr",
                                        "long_name"=>"time",
                                        "standard_name" => "time",
                                        "axis"=>"T"), 
                                            deflatelevel=DEFLATE_LEVEL) 

        defVar(ds, "step", steps, ("time",), deflatelevel=DEFLATE_LEVEL)

        # 2D Variables
        defVar(ds, "x", DTYPES["x"], ("id", "time"), 
            attrib=Dict(
                "units"         => "m", 
                "long_name"     => "x", 
                "standard_name" => "projection_x_coordinate",
                "axis"          => "X"
            ), deflatelevel=DEFLATE_LEVEL) 
        
        defVar(ds, "z", DTYPES["z"], ("id", "time"), 
            attrib=Dict(
                "units"         => "m", 
                "long_name"     => "z", 
                "standard_name" => "depth",
                "axis"          => "Z"
            ), deflatelevel=DEFLATE_LEVEL) 
            
        
        defVar(ds, "lithology", DTYPES["lithology"], ("id", "time"), 
            attrib=Dict(
                "long_name"   => "lithology",
                "coordinates" => "time z x"
            ), deflatelevel=DEFLATE_LEVEL) 
        
        # CF Global Attributes
        ds.attrib["Conventions"] = "CF-1.8"
        ds.attrib["featureType"] = "trajectory"
        
        ds.attrib["description"] = "particle trajectories"
        ds.attrib["reference_timestep"] = "$ref_step"
    end
end

# ======= processing the trajectories =======

function process_trajectories(nc_fname::String, steps::Vector{Int32}, target_ids::Vector{Int32}, id_map::Vector{Int32}, Lz::Float32, ncores::Int)
    num_particles = length(target_ids)
    num_steps = length(steps)
    
    chunk_size = ceil(Int, num_steps / CHUNKS)
    step_chunks = collect(Iterators.partition(1:num_steps, chunk_size))

    global_start = time()

    for (chunk_idx, indices) in enumerate(step_chunks)
        sub_steps = steps[indices]
        len_chunk = length(sub_steps)
        println("Running chunk $chunk_idx/$CHUNKS - Steps $(sub_steps[1]) to $(sub_steps[end])")

        
        buffer_x = fill(DTYPES["x"].(-999.0), num_particles, len_chunk)
        buffer_z = fill(DTYPES["z"].(-999.0), num_particles, len_chunk)
        buffer_litho = fill(DTYPES["lithology"].(-1), num_particles, len_chunk)
        
        chunk_start = time()
        progress_counter = Threads.Atomic{Int}(0)

        @threads for i in eachindex(sub_steps)
            step = sub_steps[i]
            
            for core in 0:(ncores-1)
                fpath = joinpath("steps", "step_$(step)_$(core).txt")
                if !isfile(fpath) continue end
                
                data = CSV.File(fpath, header=false, delim=' ', ignorerepeated=true, 
                select=[1,2,3,4], types=[DTYPES["x"], DTYPES["z"], DTYPES["id"], DTYPES["lithology"], Float64], 
                silencewarnings=true)
                
                @inbounds for row in data
                    id = row.Column3
                    
                    # Avoid OutOfBounds if a new ID appears 
                    if id >= 0 && id < length(id_map)
                        idx = id_map[id + 1]
                        
                        if idx > 0
                            buffer_x[idx, i] = row.Column1
                            buffer_z[idx, i] = row.Column2
                            buffer_litho[idx, i] = row.Column4
                        end
                    end
                end
            end
            
            Threads.atomic_add!(progress_counter, 1)
            if progress_counter[] % 10 == 0
                speed = (time() - chunk_start) / progress_counter[]
                active_ram = round(Base.gc_live_bytes() / 1000^2, digits=1) 
                @info "[particles] Progress: $(progress_counter[])/$len_chunk | Speed: $(round(speed, digits=2))s/step | Active RAM: $(active_ram)MB"
            end
        end

        Dataset(nc_fname, "a") do ds
            ds["x"][:, indices] = buffer_x
            ds["z"][:, indices] = buffer_z
            ds["lithology"][:, indices] = buffer_litho
        end
        
        buffer_x = buffer_z = buffer_litho = nothing
        GC.gc()
    end
    
    total_elapsed = time() - global_start
    @info "Finished! Total time: $(round(total_elapsed / 60, digits=2)) minutes."
end



function main()
    
    if length(ARGS) < 2
        println("Usage: julia script.jl   SCENARIO_DIR   REFERENCE_TIME_STEP")
        println("Assuming REFERENCE_TIME_STEP=0")
        ref_step = 0
    else
        ref_step = parse(Int, ARGS[end]) 
    end
    
    data_dir = ARGS[end-1]
    cd(data_dir) 

    params = read_param("param.txt")
    Lx = parse(Float32, params["lx"])
    Lz = parse(Float32, params["lz"]) 
    
    steps = get_all_steps()
    times = Float64[read_time(s) for s in steps]
    ncores = get_ncores()


    # Change this values to select particles within a box.
    global SELECT_X = (0.0f0, Lx) 
    global SELECT_Z = (-Lz, 0.0f0) 

    println("$ncores cores were used, and there are $(length(steps)) time steps.") 
    
    target_ids, id_map = extract_reference_ids(ref_step, ncores)
    
    nc_fname = "particles_trajectories.nc"
    create_particles_nc(nc_fname, target_ids, steps, times, ref_step)
    
    process_trajectories(nc_fname, steps, target_ids, id_map, Lz, ncores)
end

main()