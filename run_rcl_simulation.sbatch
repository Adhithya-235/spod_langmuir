#!/bin/bash
#SBATCH --job-name=rcl_code
#SBATCH --output=/mnt/lz01/chini/as1751-new/Research_Active_Projects/spod_langmuir/logger.out
#SBATCH --error=/mnt/lz01/chini/as1751-new/Research_Active_Projects/spod_langmuir/logger.err
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=9000

# Navigate to code directory, load module 
cd $HOME/Research_Active_Projects/spod_langmuir/
module load Anaconda3/2019.10
module load MATLAB/R2021a

# Create a folder with the current date and time
folder_name=$(date +'%Y-%m-%d_%H-%M-%S')
mkdir "$folder_name"

# Make plotting subfolders
mkdir -p "$folder_name/plots/timeseries"
mkdir -p "$folder_name/plots/frames"

# Run simulation, merge processes
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 langcirc_rclp.py --LaT=0.1 --La=0.02 --Nx=64 --Ny=64 --Nz=64 --Tend=50.00 
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs field_snapshots
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs energy_timeseries
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs timestep_tracker
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs checkpointing_data

# Move generated data and logging files to the created folder
mv field_snapshots energy_timeseries timestep_tracker checkpointing_data logger.* "$folder_name/"
cp langcirc_rclp.py run_rcl_simulation.sbatch "$folder_name/"

# Create analysis plots after navigating to post processing code folder
cd $HOME/Research_Active_Projects/spod_langmuir/post_processing/
matlab -batch "clear; close all; clc; folder_name='$folder_name'; maxs=1; plot_cwke"
matlab -batch "clear; close all; clc; folder_name='$folder_name'; svec=[1]; wrap=0; plot_fields"

# Display a message indicating completion
echo "Simulation complete. Post processing complete. Results saved in $folder_name."