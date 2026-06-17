using NCDatasets
using Glob
using CSV
using Printf
using DataFrames
using Base.Threads

# --- Funções Auxiliares (adaptadas de seus scripts) ---

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

function get_all_steps()::Vector{Int}
    pattern = joinpath("time", "time_*.txt")
    files = glob(pattern)
    
    times = [parse(Int, split(splitext(basename(f))[1], '_')[end]) for f in files]
    return sort(times)
end

function read_time(step::Int)::Float32
    time::Float32=0.0
    open(joinpath("time", "time_$step.txt")) do file
        line = readline(file)
        line = split(line,"   ")[2]
        time = parse(Float32,line)/1e6 #Myr
    return time
    end
end

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

    for idx in negative_indices
        i::Int16, j::Int16 = idx[1], idx[2]
        found = false
        for (di, dj) in [(-1, 0), (1, 0), (0, -1), (0, 1)]
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


#plotando de cabeça pra baixo
function convert_litho_to_nc(pdict::Dict{String,Any})
    nc_fname = "lithology.nc"
    Nx = pdict["Nx"]
    Nz = pdict["Nz"]
    Lx = pdict["Lx"]
    Lz = pdict["Lz"]
    steps = pdict["steps"]
    num_steps = pdict["num_steps"]
    times = pdict["times"]
    ncores = pdict["ncores"]
    litho_dict = pdict["litho_dict"]
    dfllevel::Int8 = 6 #compression level 1-9
    
    # A malha de litologia tem resolução maior
    Nxl = (Nx-1)*5 + 1
    Nzl = (Nz-1)*5 + 1
    
    full_buffer = zeros(Int8, Nxl, Nzl, num_steps)
    
    # --- Loop sobre os passos de tempo para preencher os dados ---
	#nc_lock = ReentrantLock()
	progress_counter = Threads.Atomic{Int}(0)
	total_steps = length(steps)
	start_time = time() # Start global timer
	
	#NEW IMPLEMENTS
	#n_julia_threads = Threads.nthreads()
	#thread_meshes = [fill(Int8(-1), Nzl, Nxl) for _ in 1:n_julia_threads] #each thread will be responsible for a slice
	
	@threads for i in eachindex(steps)
	    step = steps[i]
	    #@printf("Processing lithology for step: %d (%d of %d)\r", step, i, num_steps)
	    
	    # --- DEBUG: Adicione um conjunto para rastrear códigos de litologia ---
	    #unique_codes = Set{Int}()
	    # ----------------------------------------------------------------

	    # 1. Criar a malha de litologia para o passo atual
	    #litho_mesh = zeros(Int8, Nzl, Nxl) .- Int8(1)
		
		#A thread safe implementation in order to not create unnecessary litho_mesh's in memory 
		#tid = Threads.threadid()
		#litho_mesh = thread_meshes[tid]
		litho_mesh = fill(Int8(-1), Nzl, Nxl)
		#fill!(litho_mesh, Int8(-1)) #reset with -1
		
	    for core in 0:(ncores-1)
		fpath=joinpath("lithos","litho_$(step)_$core.txt")
			if isfile(fpath)
				x_core, z_core, litho_core = read_litho_file(fpath)
				for k in eachindex(x_core)
				    litho_mesh[z_core[k] + 1, x_core[k] + 1] = litho_core[k]
				end
			end
	    end

	    # 2. Processar a malha (preencher vazios e mapear valores)
	    replace_negatives_with_neighbors!(litho_mesh)
	    litho_mesh .= get.(Ref(litho_dict), litho_mesh, litho_mesh)

	    # Store finished slice into our 3D buffer
	    full_buffer[:, :, i] = reverse(litho_mesh, dims=1)'
	    Threads.atomic_add!(progress_counter, 1)
	    if progress_counter[] % 10 == 0
        @info "Progress: $(progress_counter[]) / $total_steps ($(round(progress_counter[]/total_steps*100, digits=1))%)"
    	end
	end
    total_elapsed = time() - start_time
    @info "Finished! Total time: $(round(total_elapsed / 60, digits=2)) minutes."
    
    x_coords_litho = range(0.0f0, Lx, length=Nxl)
    z_coords_litho = range(0.0f0, Lz, length=Nzl)
    
    
    Dataset(nc_fname, "c") do ds
        # --- Definir Dimensões e Variáveis ---
        defDim(ds, "time", num_steps)
        defDim(ds, "x", Nxl)
        defDim(ds, "z", Nzl)

        defVar(ds, "time", Float32.(times), ("time",), 
            attrib=Dict("units"=>"Myr", "long_name"=>"Model time","axis"=>"T"),
            deflatelevel=dfllevel, shuffle=true)
        defVar(ds, "x", x_coords_litho, ("x",), 
            attrib=Dict("units"=>"m", "long_name"=>"x-coordinate", "axis"=>"X"),
            deflatelevel=dfllevel, shuffle=true)
        defVar(ds, "z", z_coords_litho, ("z",), 
            attrib=Dict("units"=>"m", "long_name"=>"z-coordinate", "axis"=>"Z"),
            deflatelevel=dfllevel, shuffle=true)
        
        defVar(ds, "lithology", Int8, ("x", "z", "time"),
            attrib=Dict(
                "long_name"=>"Lithology",
                "comment"=>"Codes are mapped: 1=Mantle, 2=LC, 3=UC, 4=Salt, SEDIMENTS, 12: Air"
            ),
            deflatelevel=dfllevel, shuffle=true)

        ds["lithology"][:, :, :] = full_buffer
        println("\nSaved to $nc_fname")
    end
end

data_dir = ARGS[1]
cd(data_dir)

params = read_param("param.txt")
Nx = parse(Int, params["nx"])
Nz = parse(Int, params["nz"])
Lx = parse(Float32, params["lx"])
Lz = parse(Float32, params["lz"])

steps::Vector{Int} = get_all_steps()
num_steps::Int = length(steps)
times = Float32[read_time(step) for step in steps]
#name::String = split(basename(glob(joinpath("lithos","litho_0_*.txt"))),'.')[1]
ncores::Int = size(glob(joinpath("lithos","litho_0_*.txt")))[1]
println("Encontrados $num_steps passos de tempo, de $(steps[1]) a $(steps[end]).")
println("Número de núcleos detectados: $ncores")

#mudar litologia com base nas camadas do interfaces
litho_dict = Dict{Int8,Int8}(
    2 => 1,  # Lithospheric mantle
    3 => 1,  # Lithospheric mantle
    4 => 2,  # Lower crust
    5 => 3,  # Upper crust
    6 => 4,  # Salt
    7 => 5,  # Sed 0
    8 => 6,  # Sed 1
    9 => 7,
    10 => 8,
    11 => 9,
    12 => 10,
    13 => 11,  #sed pre-rift
    14 => 12  #air
    )


pdict = Dict{String,Any}(
    "Nx"=>Nx, "Nz"=>Nz, "Lx"=>Lx, "Lz"=>Lz,
    "steps"=>steps, "num_steps"=>num_steps, "times"=>times,
    "ncores"=>ncores, "litho_dict"=>litho_dict
)

println("Convertendo dados de litologia...")
convert_litho_to_nc(pdict)
println("\nConversão concluída!")
