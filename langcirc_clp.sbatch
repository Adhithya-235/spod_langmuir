#!/bin/bash
#SBATCH --job-name=la2e-2
#SBATCH --output=/mnt/lz01/chini/as1751-new/Research_Active_Projects/Craik_Leibovich_DNS/codebcfix/la%j.out
#SBATCH --error=/mnt/lz01/chini/as1751-new/Research_Active_Projects/Craik_Leibovich_DNS/codebcfix/la%j.err
#SBATCH --mail-user=as1751@wildcats.unh.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=4
#SBATCH --time=9000

cd $HOME/Research_Active_Projects/Craik_Leibovich_DNS/codebcfix/
module load Anaconda3/2019.10
mpiexec -n 96 $HOME/.conda/envs/dedalus/bin/python3 langcirc_clp.py --LaT=0.1 --Re=100.0 --Del=1.0 --Nx=128 --Ny=128 --Nz=64 --Tend=600.00 --La1=1 --La2=-4
mpiexec -n 96 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs la1e-4
mpiexec -n 96 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs la1e-4a
mpiexec -n 96 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs la1e-4t
mpiexec -n 96 $HOME/.conda/envs/dedalus/bin/python3 -m dedalus merge_procs la1e-4p
