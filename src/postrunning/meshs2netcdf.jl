using NCDatasets
using Glob
using CSV
using Printf
using DataFrames
using Base.Threads
using StatsBase

global AIR_DENSITY_THRESHOLD = 20
global LITHOLOGY_DATATYPE = UInt8
global VARIABLES = ["density", "viscosity", "pressure", "strain","strain_rate","temperature","velocity","surface","heat"]

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
    "temperature"=>"°C",
    "velocity"=>"m/s",
    "surface"=>"m",
    "heat"=>"W/kg",
    "thermal_diffusivity"=>"m2/s",
    "Phi"=>"dimensionless",
    "dPhi"=>"1/s",  # Need to confirm
    "X_depletion"=>"dimensionless"
)

global DTYPES =Dict{String,DataType}(
    "x"=>Float32,
    "y"=>Float32,
    "z"=>Float32,
    "time"=>Float32,

    "density"=>UInt32,
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
    "X_depletion"=>Float64
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
    Ly::Int
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

function read_data(var::String, step::Int, nxnz::Tuple; veloc::Bool=false, surface::Bool=false,vartype::DataType=Float64)
    
    file = "$(var)/$(var)_$(step).txt"
    skipto::Int = 0
    if surface skipto=0 else skipto=3 end
    C = CSV.File(file, header=false, comment="P", skipto=skipto,types=vartype)|>CSV.Tables.matrix
    
    Nx,Nz = nxnz
    if veloc 
        C[C .< 1e-200] .= 0
        vx = transpose(reshape(C[1:2:end], (Nx, Nz)))
        vy = transpose(reshape(C[2:2:end], (Nx, Nz)))
        R = (vx,vy)
        return R
    
    elseif surface
        return (C[1:end,1],C[1:end,2]) #sx, sy
    
    else
        C[C .< 1e-200] .= 0
        R = transpose(reshape(C[1:end], (Nx, Nz)))
        return R
    end
end

function get_all_steps()::Vector{Int}
    pattern = joinpath("time", "time_*.txt")
    files = glob(pattern)
    
    steps = Int32[parse(Int32, split(splitext(basename(f))[1], '_')[end]) for f in files]
    return sort(steps)
end

function read_time(step::Int)::Float64
    time::Float64=0.0
    open(joinpath("time", "time_$step.txt")) do file
        line = readline(file)
        line = split(line,"   ")[2]
        time = parse(Float64,line)/1e6 #Myr
    return time
    end
end

function converter(variable::String, scen::MandyocScenario, mesh::mesh2D)
    nc_fname = "$variable.nc"
    Nx = mesh.Nx
    Nz = mesh.Nz
    Lx = mesh.Lx
    Lz = mesh.Lz

    dtypes = scen.datatypes[variable]
    vtype = dtypes[variable]
    x_coords = dtypes["x"].(range(0.0f0, Lx, length=Nx))
    z_coords = dtypes["z"].(range(0.0f0, Lz, length=Nz))

    steps = scen.steps
    times = scen.times
    num_steps = length(steps)

    units = scen.units
    
    dfllevel::Int8 = 7 # compression level 1-9


    veloc = (variable == "velocity")
    surface = (variable == "surface")
    
    local buffer_vx, buffer_vy, buffer_var, buffer_surf, surface_nx, surface_x_coords
    
    if surface
        sx_sample,_ = read_data("surface", 0, (Nx,Nz), veloc=false, surface=true)
        surface_nx = size(sx_sample)[1]
        surface_x_coords = dtypes["x"].( range(0.0f0, Lx, length=surface_nx) )
		buffer_surf = zeros(dtypes["surface"], surface_nx, num_steps)

    elseif veloc
        buffer_vx = zeros(dtypes["velocity"], Nx, Nz, num_steps)
        buffer_vy = zeros(dtypes["velocity"], Nx, Nz, num_steps)
    else
        buffer_var = zeros(vtype, Nx, Nz, num_steps)
    end
    
    #Multithreading processing
    progress_counter = Threads.Atomic{Int}(0)
    start_time = time()
    
    @threads for i in eachindex(steps)
        step = steps[i]

        data = read_data(variable,step,(Nx,Nz),veloc=veloc,surface=surface,vartype=vtype)
        if data !== nothing
            
	    if veloc
                dens = read_data("density",step,(Nx,Nz),veloc=false, surface=false)
                vx,vy = data
                vx[dens.<AIR_DENSITY_THRESHOLD] .= 0
                vy[dens.<AIR_DENSITY_THRESHOLD] .= 0
		        buffer_vx[:, :, i] = vx'
                buffer_vy[:, :, i] = vy'
                
            elseif surface

                sx,sy = data

                buffer_surf[:, i] = sy
            else
                dens = read_data("density",step,(Nx,Nz),veloc=false, surface=false)
                data[dens.<AIR_DENSITY_THRESHOLD] .= 0
                buffer_var[:, :, i] = data'
            end
	    
            else
                @warn "No data found for $fpath at step $step"
            end
        
        #Tracker
			Threads.atomic_add!(progress_counter, 1)
			if progress_counter[] % 10 == 0
				speed = (time() - start_time) / progress_counter[]
				@info "[$variable] Progress: $(progress_counter[])/$num_steps | Speed: $(round(speed, digits=2))s/step"
			end
        
        end
    
    Dataset(nc_fname,"c") do ds #criar o arquivo nc

        defDim(ds,"time",num_steps)
        defVar(ds,"time", times, ("time",),
        attrib=Dict("units"=>units["time"],"long_name"=>"Time","axis"=>"T"),
                                                                deflatelevel=dfllevel, shuffle=true)
        
        defVar(ds,"steps", steps, ("time",),
        attrib=Dict("units"=>"","long_name"=>"Time steps"),
                deflatelevel=dfllevel, shuffle=true)
                                                         
        if veloc
            defDim(ds,"x",Nx)
            defDim(ds,"z",Nz)

            defVar(ds,"x",x_coords,("x",),attrib=Dict("units"=>units["x"],"long_name"=>"x","axis" => "X"),
                                                                deflatelevel=dfllevel, shuffle=true,
                                                                )
            defVar(ds,"z",z_coords,("z",),attrib=Dict("units"=>units["z"],"long_name"=>"z","axis"=>"Z"),
                                                                deflatelevel=dfllevel, shuffle=true)
            
            defVar(ds,"vx",vtype,("x","z","time"),attrib=Dict("units"=>units[variable],
                                                                "long_name"=>"vx"),
                                                                deflatelevel=dfllevel, shuffle=true)
            
            defVar(ds,"vy",vtype,("x","z","time"),attrib=Dict("units"=>units[variable],
                                                                "long_name"=>"vy",),
                                                               deflatelevel=dfllevel, shuffle=true)
            ds["vx"][:, :, :] = buffer_vx
            ds["vy"][:, :, :] = buffer_vy
            
        elseif surface
            defDim(ds,"x",surface_nx)
            defVar(ds,"x",surface_x_coords,("x",),attrib=Dict("units"=>units["x"],"long_name"=>"x","axis" => "X"), 
            													deflatelevel=dfllevel, shuffle=true)

            defVar(ds,variable,vtype,("x","time"),attrib=Dict("units"=>units[variable],"long_name"=>variable),
                                                                deflatelevel=dfllevel, shuffle=true)
            ds[variable][:, :] = buffer_surf
            
        else
            defDim(ds,"x",Nx)
            defDim(ds,"z",Nz)
            defVar(ds,"x",x_coords,("x",),attrib=Dict("units"=>units["x"],"long_name"=>"x-coordinate","axis" => "X"), 
            													deflatelevel=dfllevel,shuffle=true)
            defVar(ds,"z",z_coords,("z",),attrib=Dict("units"=>units["z"],"long_name"=>"z-coordinate","axis"=>"Z"),
                                                                deflatelevel=dfllevel,shuffle=true)
            defVar(ds,variable,vtype,("x","z","time"),attrib=Dict("units"=>units[variable],"long_name"=>variable),
                                                                deflatelevel=dfllevel, shuffle=true)
                                                                
            ds[variable][:, :, :] = buffer_var
        end
	
        
        println("\nSaved to $nc_fname")
    end
end

function build_scenario(params::Dict)

    dims = parse(Int, get(params, "dimensions", "2")) # if the model is 2d or 3d

    # Change data types in case of non dimensional scenarios 
    if get(params,"iterative","direct") == "iterative" || get(params,"nondimensionalization","False") == "True" || dims == 3
        for v in keys(DTYPES) DTYPES[v] = Float64 end
        global AIR_DENSITY_THRESHOLD = 0
        println("You did run a non dimensional model, all datatypes changed to Float64.")
    end

    Nx = parse(Int, params["nx"]) # elements in x
    Nz = parse(Int, params["nz"]) # elements in z
    Lx = parse(DTYPES["x"], params["lx"]) # x length
    Lz = parse(DTYPES["z"], params["lz"]) # z length
    thick_air::DTYPES["z"] = 40e3 # m

    if dims == 3 
        Ny = parse(Int, params["ny"])
        Ly = parse(Int, params["ly"])
        mesh = mesh3D(Nx,Ny,Nz,Lx,Ly,Lz)
    else
        mesh = mesh2D(Nx,Nz,Lx,Lx)
    end

    # Finding all steps
    steps = get_all_steps()
    times = DTYPES["time"][read_time(step) for step in steps]
    CSV.write("times.csv", DataFrame(step=steps, time_myr=times))
    println("$(length(steps)) time steps were found, from $(times[1]) Myr [$(steps[1])] to $(times[end]) Myr [$(steps[end])].")

    return MandyocScenario(dims, steps, times, thick_air, UNITS, DTYPES), mesh
end

function main()
    
data_dir = ARGS[end] # Scenario directory
cd(data_dir)

# Basic parameters
params = read_param("param.txt")

if get(params, "magmatism", "off") == "on"
    push!(VARIABLES, "Phi")
    push!(VARIABLES, "dPhi")
    push!(VARIABLES, "X_depletion")
    println("magmatism=on")
end

if get(params, "export_thermal_diffusivity", "False") == "True"
    push!(VARIABLES, "thermal_diffusivity")
end

scen, mesh = build_scenario(params)

for var in unique(VARIABLES)
    println("Converting: $var")
    converter(var,scen,mesh)
end

println("All variables were converted to NetCDF4")
println("Finished")

end

main()