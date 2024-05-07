#!/bin/bash
#SBATCH --job-name=la2e-2
#SBATCH --output=/mnt/lz01/chini/as1751-new/Research_Active_Projects/spod_langmuir/la%j.out
#SBATCH --error=/mnt/lz01/chini/as1751-new/Research_Active_Projects/spod_langmuir/la%j.err
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=9000

cd $HOME/Research_Active_Projects/Craik_Leibovich_DNS/codebcfix/
module load Anaconda3/2019.10
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 langcirc_clp.py --LaT=0.1 --La=0.02 --Nx=64 --Ny=64 --Nz=64 --Tend=600.00 --La1=2 --La2=-2
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs la2e-2
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs la2e-2a
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs la2e-2t
mpiexec -n 16 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs la2e-2p
