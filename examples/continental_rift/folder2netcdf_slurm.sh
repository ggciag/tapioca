#!/usr/bin/env bash

#SBATCH --job-name=TAPIOCA_NC  
#SBATCH --output=slurm-%A_%a.log   
#SBATCH --mail-user=jbueno@usp.br
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=1                 
#SBATCH --cpus-per-task=16         
#SBATCH --time=72:00:00            
#SBATCH --mem=20G

# --- Settings ---
current_scenario=$1
tapioca_dir="$(PATH_TO_TAPIOCA)/tapioca/src/postrunning"
cenario_dir=$1


# --- Executing ---
echo "========================================================"
echo "Processing Scenario: $cenario_dir"
echo "Running on: $(hostname) with $SLURM_CPUS_PER_TASK CPUs"
echo "========================================================"

cd "$cenario_dir" || exit 1
bash "$tapioca_dir/0_organize_outputs.sh"
julia -t $SLURM_CPUS_PER_TASK "$tapioca_dir/meshes2netcdf.jl" "$cenario_dir"
julia -t $SLURM_CPUS_PER_TASK "$tapioca_dir/particles2netcdf.jl" "$cenario_dir" 0

echo "NETCDF conversion is concluded: $current_scenario. Zipping and removing dirs"

data_dirs="heat viscosity velocity density pressure temperature surface strain_rate strain thermal_diffusivity time"

echo "Zstd compression"

rm -rf "$cenario_dir"/lithos

tar --use-compress-program="zstd -8 -T$SLURM_CPUS_PER_TASK" \
    -cf "$cenario_dir/raw_data.tar.zst" \
    --remove-files \
    ${data_dirs}

tar --use-compress-program="zstd -8 -T$SLURM_CPUS_PER_TASK" \
    -cf "$cenario_dir/raw_steps.tar.zst" \
    --remove-files \
    "steps"

echo "Processing of $current_scenario finished."
