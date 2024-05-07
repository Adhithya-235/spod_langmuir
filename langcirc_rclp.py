#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage: 
  langcirc_clp.py [--LaT=<turbulent_langmuir> --La=<laminar_langmuir> \
  --Nx=<Nx> --Ny=<Ny> --Nz=<Nz> --Tend=<stop_time>] 
  
Options:
  --LaT=<turbulent_langmuir>        Turbulent Langmuir Number [default: 0.1]
  --La=<laminar_langmuir>           Friction Reynolds Number [default: 0.02]
  --Nx=<Nx>                         Number of downwind modes [default: 64]
  --Ny=<Ny>                         Number of crosswind modes [default: 64]
  --Nz=<Nz>                         Number of vertical modes [default: 64]
  --Tend=<stop_time>                Simulation stop time [default: 400.0]  
"""

"""
Dedalus script for the direct numerical simulation of the reduced Craik-Leibovich 
equations as presented in Chini et al (2009).

This script uses a Fourier basis in x and y (PBCs) and a Chebyshev basis in z 
with nonhomogeneous Neumann BCs for the horizontal velocity 'u', homogeneous 
Neumann BCs for the cross-wind velocity 'v' and Dirichlet BCs for the vertical
velocity 'w'. 

This script can be run serially or in parallel, and uses the built-in analysis 
framework to save data in HDF5 files. The 'merge_procs' command can be used to 
merge distributed analysis sets from parallel runs. 

To run and merge using 4 processes, for instance, you could use: 
    $ mpiexec -n 4 python3 langcirc_clp.py
    $ mpiexec -n 4 python3 -m dedalus merge_procs rb1
    
This script can restart the simulation from the last save of the original
output to extend the integration. This requires that the output files from the
original simulation are merged, and the last is symlinked or copied to 
'restart.h5'

"""
# SET UP ENVIRONMENT

import pathlib
import time
import h5py
import numpy as np
from mpi4py import MPI
from docopt import docopt
from dedalus import public as de
from dedalus.extras import flow_tools
import logging
logger = logging.getLogger(__name__)

# DEFINE PARAMETERS

args = docopt(__doc__)
LangmuirT  = float(args['--LaT'])                                                             # Turbulent Langmuir Number
Langmuir   = float(args['--La'])                                                              # Wave Reynolds Number
Lx, Ly, Lz = (4.0*np.pi, np.pi, 1.0)                                                          # Box Size
Nx, Ny, Nz = (int(args['--Nx']), int(args['--Ny']), int(args['--Nz']))                        # No. of Gridpoints
stop_time  = float(args['--Tend'])                                                            # Sim. stop time
la1        = int(args['--La1'])                                                               # File naming sequence, first number
la2        = int(args['--La2'])                                                               # File naming sequence, second number

# LOGGER: RECORD INPUT PARAMETERS

if MPI.COMM_WORLD.rank == 0:
    logger.info("Running rCL simulation for LaT={:.3e}, La={:.3e}".format(LangmuirT, Langmuir))
    
# CREATE BASES AND DOMAIN

x_basis = de.Fourier("x", Nx, interval=(0, Lx), dealias=3/2)
y_basis = de.Fourier("y", Ny, interval=(0, Ly), dealias=3/2)
z_basis = de.Chebyshev("z", Nz, interval=(-1.0*Lz, 0), dealias=3/2)
domain  = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64, mesh=[4,4])
x, y, z    = domain.grids(scales=1)
xd, yd, zd = domain.grids(scales=3/2)

# PROBLEM SETUP

problem = de.IVP(domain, variables=['U', 'V', 'W', 'P', 'Uz', 'Vz', 'Wz'])
problem.meta['Uz', 'Vz', 'W', 'P']['z']['dirichlet'] = True

# NON-CONSTANT COEFFICIENTS (STOKES DRIFT VELOCITY)

ncc = domain.new_field(name='Us')
ncc['g'] = 1.0 + z
ncc.meta['x']['constant'] = True
ncc.meta['y']['constant'] = True

# EQUATION ENTRY SUBSTITUTIONS

problem.parameters['La']   = Langmuir
problem.parameters['LaT']  = LangmuirT
problem.parameters['Us']   = ncc
problem.parameters['Lx']   = Lx
problem.parameters['Ly']   = Ly
problem.parameters['Lz']   = Lz

problem.substitutions['Lap(A, Az)']       = 'dy(dy(A)) + dz(Az)'
problem.substitutions['Adv(A, B, C, Cz)'] = 'A*dy(C) + B*Cz'

# ANALYSIS SUBSTITUTIONS

#---- OPERATORS

problem.substitutions['xavg(A)'] = 'integ(A, "x")/Lx'
problem.substitutions['yavg(A)'] = 'integ(A, "y")/Ly'
problem.substitutions['zavg(A)'] = 'integ(A, "z")/Lz'
problem.substitutions['Vavg(A)'] = 'integ(integ(integ(A, "x"), "y"), "z")/(Lx*Ly*Lz)'

#---- VORTICITY

problem.substitutions['omega_x'] = 'dy(W) - Vz'

#---- ENERGIES

problem.substitutions['KE']   = 'Vavg(((LaT**(8.0/3.0))*(La**(-2.0/3.0))*U*U + V*V + W*W)/2.0)'
problem.substitutions['UNKE'] = 'Vavg((U*U + V*V + W*W)/2.0)'
problem.substitutions['CWKE'] = 'Vavg((V*V + W*W)/2.0)'

# EVOLUTION EQUATIONS 
    
problem.add_equation("dt(U) + dx(P) - La*Lap(U, Uz) = -Adv(V, W, U, Uz)") 
problem.add_equation("dt(V) + dy(P) - La*Lap(V, Vz) - Us*(dy(U) - dx(V)) = - Adv(V, W, V, Vz)")
problem.add_equation("dt(W) + dz(P) - La*Lap(W, Wz) - Us*(Uz - dx(W)) = - Adv(V, W, W, Wz)")

# CONSTRAINT EQUATIONS
    
problem.add_equation("dy(V) + Wz = 0")
problem.add_equation("dz(F) = 0", condition="(ny!=0)")
problem.add_equation("P - dz(F) = 0", condition="(ny==0)")

# AUXILIARY EQUATIONS (DEFINE z-DERIVATIVES)

problem.add_equation("Uz - dz(U) = 0")
problem.add_equation("Vz - dz(V) = 0")
problem.add_equation("Wz - dz(W) = 0")

# BOUNDARY CONDITIONS

problem.add_bc("left(Uz) = 1")
problem.add_bc("right(Uz) = 1")
problem.add_bc("left(Vz) = 0")
problem.add_bc("right(Vz) = 0")
problem.add_bc("left(W) = 0")
problem.add_bc("right(W) = 0", condition="(ny != 0)")
problem.add_bc("left(F) = 0")
problem.add_bc("right(F) = 0", condition="(ny == 0)")

# BUILD SOLVER

ts = de.timesteppers.SBDF3
solver = problem.build_solver(ts)
logger.info('Solver built')

# INITIAL CONDITIONS OR RESTART CODE

if not pathlib.Path('restart.h5').exists():
    
    U   = solver.state['U']
    Uz  = solver.state['Uz']
    
    # RANDOM PERTURBATIONS, INITIALIZED GLOBALLY FOR SAME RESULTS IN PARALLEL
    
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=23)
    noise = rand.standard_normal(gshape)[slices]
    
    # DAMP PERTURBATIONS AT WALLS
    
    zb, zt = z_basis.interval
    pert = 1e-3*noise*(zt-z)*(z-zb)
    
    # INITIALIZE Uz, INTEGRATE TO GET U
    
    Uz['g'] = 1.0 + pert
    Uz.antidifferentiate('z', ('left', 0), out=U)

    # INTEGRATION PARAMETERS
    
    dt = 1e-5
    solver.stop_sim_time = stop_time
    fh_mode = 'overwrite'
    
else:
    
    # RESTART

    write, last_dt = solver.load_state('restart.h5', -1)
    
    # INTEGRATION PARAMETERS

    dt = last_dt
    solver.stop_sim_time = 700
    fh_mode = 'append'
    
# ANALYSIS

snapshot = solver.evaluator.add_file_handler("field_snapshots", sim_dt=0.2, max_writes=400, mode=fh_mode)
snapshot.add_task("U", name = 'U')
snapshot.add_task("V", name = 'V')
snapshot.add_task("W", name = 'W')
snapshot.add_task("P", name = 'P')
snapshot.add_task("omega_x", name = 'O')

globalp = solver.evaluator.add_file_handler("energy_timeseries", sim_dt=0.02, max_writes=10000000, mode=fh_mode)
globalp.add_task("KE", name = 'KE')
globalp.add_task("CWKE", name = 'CWKE')

globalt = solver.evaluator.add_file_handler("timestep_tracker", iter=1, max_writes=10000000, mode=fh_mode)
globalt.add_task("CWKE", name = 'CWKE')

photo = solver.evaluator.add_file_handler("checkpointing_data", sim_dt=0.4, max_writes=100, mode=fh_mode)
photo.add_system(solver.state)
photo.add_task("Uz", name = 'Uz')

# CFL 

CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.3, max_change=1.5, min_change=0.5, max_dt=0.01, threshold=0.1)
CFL.add_velocity('V', 1)
CFL.add_velocity('W', 2)

# FLOW TOOLS

flow = flow_tools.GlobalFlowProperty(solver, cadence = 100)
flow.add_property("dy(V) + Wz", name="divUP")
flow.add_property("KE", name="TKE")

# MAIN LOOP

try:  
    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        
        dt = CFL.compute_dt()
        solver.step(dt)
        if (solver.iteration-1) % 100 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e, Div(V,W): %e, ToKE: %e' %(solver.iteration, solver.sim_time, dt, flow.max('divUP'), flow.max('TKE')))     
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))


