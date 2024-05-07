#!/bin/bash

# Create a folder with the current date and time
folder_name = $(date + '%Y-%m-%d_%H-%M-%S')
mkdir "$folder_name"

# Send the simulation to a remote cluster, merge processes
sbatch langcirc_rclp.sbatch
echo "Job sent to cluster."

# Move generated data and logging files to the created folder
mv field_snapshots energy_timeseries timestep_tracker checkpointing_data langcirc_rclp.py langcirc_rclp.sbatch "$folder_name/"

# Display a message indicating completion
echo "Simulation complete. Results saved in $folder_name."
