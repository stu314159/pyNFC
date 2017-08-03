#!/p/home/sblair/anaconda2/bin/python

#pyNFC_preprocess.py
"""
   pre-processing steps for pyNFC code.  Basic geometry definition
   and partitioning information created and input file created prior to
   parallel simulation run

"""

#import FluidChannel as fc
#import pyPartition as pp
#from pymetis import part_graph #<-- requires that the PrgEnv-intel module be selected
#import numpy as np
import scipy.io
import math
import argparse

parser = argparse.ArgumentParser(prog='pyNFC_preprocess.py',
                                 description='pre-processing script for pyNFC')

parser.add_argument('geom_filename',type=str)
parser.add_argument('lattice_type',type=str)
parser.add_argument('dynamics',type=int)
parser.add_argument('partition_style',type=str)
parser.add_argument('numProcs',type=int)
parser.add_argument('Num_ts',type=int)
parser.add_argument('ts_rep_freq',type=int)
parser.add_argument('Warmup_ts',type=int)
parser.add_argument('plot_freq',type=int)
parser.add_argument('Re',type=float)
parser.add_argument('dt',type=float)
parser.add_argument('Cs',type=float)
parser.add_argument('Restart_flag',type=int)
# parse input arguments
args = parser.parse_args()

# assign to required variables
geom_filename = args.geom_filename
lattice_type = args.lattice_type
dynamics = args.dynamics
partition_style = args.partition_style
numProcs = args.numProcs


#lattice_type = 'D3Q15' # [ 'D3Q15' | 'D3Q19' | 'D3Q27' ]
#partition_style = 'metis' # [ '1D' | '3D' | 'metis']

Num_ts = args.Num_ts
ts_rep_freq = args.ts_rep_freq
Warmup_ts = args.Warmup_ts
plot_freq = args.plot_freq
Re = args.Re
dt = args.dt
Cs = args.Cs
Restart_flag = args.Restart_flag

#numProcs = 6  #<--- for this version, I will need to know how many partitions I intend to create

# --- do input file processing as with genInput.py - will also add in the partitioning information ---- 
# ---- this means, I will need to know the number of processes in advance, I guess ----- though
# ---- this could alternatively be done from the pyNFC code...think about that...
#----You should not have to edit anything below this line -------------------

geom_input = scipy.io.loadmat(geom_filename)
# overall domain dimensions
Lx_p = float(geom_input['Lx_p'])
Ly_p = float(geom_input['Ly_p'])
Lz_p = float(geom_input['Lz_p'])
Lo = float(geom_input['Lo'])
Ny_divs = int(geom_input['Ny_divs'])
rho_p = float(geom_input['rho_p'])
nu_p = float(geom_input['nu_p'])
pRef_idx = int(geom_input['pRef_idx'])
#snl = list((geom_input['snl']).flatten())
#inl = list((geom_input['inl']).flatten()) # must be inlet on Z-min
#onl = list((geom_input['onl']).flatten()) # must be outlet on Z-max

Ny = math.ceil((Ny_divs-1)*(Ly_p/Lo))+1
Nx = math.ceil((Ny_divs-1)*(Lx_p/Lo))+1
Nz = math.ceil((Ny_divs-1)*(Lz_p/Lo))+1
nnodes = Nx*Ny*Nz


# non-dimensionalization
Uo = nu_p*Re/Lo
To = Lo/Uo
Uavg = Uo

Ld = 1.; Td = 1.; Ud = (To/Lo)*Uavg;
nu_d = 1./Re
dx = 1./(Ny_divs - 1.)
u_lbm = (dt/dx)*Ud
nu_lbm = (dt/(dx**2))*nu_d
omega = 1./(3.*nu_lbm+0.5)

u_conv_fact = (dt/dx)*(To/Lo)
t_conv_fact = (dt*To)
l_conv_fact = dx*Lo
#p_conv_fact = (((l_conv_fact/t_conv_fact)**2)*(1./3.))*(l_conv_fact**3)
p_conv_fact = (((l_conv_fact/t_conv_fact)**2)*(1./3.))

rho_lbm = rho_p*(l_conv_fact**3)

print 'l_conv_fact = %g.\n'%l_conv_fact
print 'p_conv_fact = %g.\n'%p_conv_fact


print 'Number of lattice points = %d.' % nnodes
print 'Number of time-steps = %d.' % Num_ts
print 'LBM viscosity = %g.' % nu_lbm
print 'LBM relaxation parameter (omega) = %g.' % omega
print 'LBM flow Mach number = %g. ' % u_lbm
print 'Nx = %d' % Nx
print 'Ny = %d' % Ny
print 'Nz = %d' % Nz

#run_dec = raw_input('Would you like to continue? [Y/n]: ')
run_dec = 'y' # just let it run

if run_dec!='n' and run_dec!='N':
    print 'Ok! Cross your fingers!!'
    # write the input file
    params = open('params.lbm','w')
    params.write('%s \n'% lattice_type) # lattice selection (keep.  We might actually use this)
    params.write('%d \n'% dynamics)
    params.write('%d \n'%Num_ts)
    params.write('%d \n'%ts_rep_freq)
    params.write('%d \n'%Warmup_ts)
    params.write('%d \n'%plot_freq)
    params.write('%g \n'%Cs)
    params.write('%g \n'%rho_lbm) # density
    params.write('%g \n'%u_lbm) # scaled maximum velocity
    params.write('%g \n'%omega) # relaxation parameter
    params.write('%d \n'%Nx) # number of nodes in the x, y and z direction
    params.write('%d \n'%Ny)
    params.write('%d \n'%Nz)
    params.write('%d \n'%Restart_flag) # 1 = load restart data; 0 = no restart
    
    # the following will not be used by the MPI code, but will be available
    # for the post-processing script
    
    params.write('%f \n'%Lx_p) # physical dimensions in the x,y and z dimensions
    params.write('%f \n'%Ly_p)
    params.write('%f \n'%Lz_p)
    
    params.write('%15.14f \n'%t_conv_fact)  # time, length and pressure conversion factors
    params.write('%15.14f \n'%l_conv_fact)
    params.write('%g \n'%p_conv_fact)
    params.write('%d \n'%pRef_idx);
    
    params.close()
    
else:
    print 'Run aborted.  Better luck next time!'



