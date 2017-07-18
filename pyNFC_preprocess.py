#pyNFC_preprocess.py
"""
   pre-processing steps for pyNFC code.  Basic geometry definition
   and partitioning information created and input file created prior to
   parallel simulation run

"""

import FluidChannel as fc
import pyPartition as pp
from pymetis import part_graph #<-- requires that the PrgEnv-intel module be selected
import numpy as np
import scipy
import math
import argparse

parser = argparse.ArgumentParser(prog='pyNFC_preprocess.py',
                                 description='pre-processing script for pyNFC')

parser.add_argument('lattice_type',type=str)
parser.add_argument('partition_style',type=str)
parser.add_argument('numProcs',type=int)
# parse input arguments
args = parser.parse_args()

# assign to required variables
lattice_type = args.lattice_type
partition_style = args.partition_style
numProcs = args.numProcs

# generate a suitable geometry:

#overall channel dimensions
aLx_p = 6.4
aLy_p = 3.0
aLz_p = 14.0
aNdivs = 11

# wall mounted brick parameters
x_c = 3.5;
z_c = 3.2;
W = 1.;
H = W;
L = W;

myObst = fc.WallMountedBrick(x_c,z_c,L,W,H);

myChan = fc.FluidChannel(Lx_p=aLx_p,Ly_p=aLy_p,Lz_p=aLz_p,obst=myObst,
                         N_divs=aNdivs)                         

# write the mat file
geom_file_stub='wall_mounted_brick';
myChan.write_mat_file(geom_file_stub);



# set simulation parameters (as used in genInput.py):
geom_filename = geom_file_stub + '.mat'
#lattice_type = 'D3Q15' # [ 'D3Q15' | 'D3Q19' | 'D3Q27' ]
#partition_style = 'metis' # [ '1D' | '3D' | 'metis']
Num_ts = 20001
ts_rep_freq = 50
Warmup_ts = 0
plot_freq = 2000
Re = 5
dt = 0.005
Cs = 0
Restart_flag = 0

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
snl = list((geom_input['snl']).flatten())
inl = list((geom_input['inl']).flatten()) # must be inlet on Z-min
onl = list((geom_input['onl']).flatten()) # must be outlet on Z-max

Ny = math.ceil((Ny_divs-1)*(Ly_p/Lo))+1
Nx = math.ceil((Ny_divs-1)*(Lx_p/Lo))+1
Nz = math.ceil((Ny_divs-1)*(Lz_p/Lo))+1
nnodes = Nx*Ny*Nz

# compute geometric data only once
x = np.linspace(0.,Lx_p,Nx).astype(np.float32);
y = np.linspace(0.,Ly_p,Ny).astype(np.float32);
z = np.linspace(0.,Lz_p,Nz).astype(np.float32);
numEl = Nx*Ny*Nz
Y,Z,X = np.meshgrid(y,z,x);

XX = np.reshape(X,int(numEl))
YY = np.reshape(Y,int(numEl))
ZZ = np.reshape(Z,int(numEl))

print 'There are %d nodes in the solid node list'%len(snl)
print 'Writing those nodes to file'
# now write this obstList to file.
obstFilename = 'snl.lbm'
obstFile = open(obstFilename,'w')
obstFile.write('%i \n'%len(snl))
for i in range(len(snl)):
    nd = int(snl[i]); nd=nd;# make snl node numbers 0-based
    obstFile.write('%i \n'% nd) 
obstFile.close()

print 'There are %d nodes in the inlet node list'%len(inl)
print 'Writing those nodes to file'
inletFileName = 'inl.lbm'
inletFile = open(inletFileName,'w')
inletFile.write('%i \n'%len(inl))
for i in range(len(inl)):
    nd = int(inl[i]); nd = nd;#make inl node numbers 0-based
    inletFile.write('%i \n'% nd) 
inletFile.close()

print 'There are %d nodes in the outlet node list'%len(onl)
print 'Writing those nodes to file'
outletFileName = 'onl.lbm'
outletFile = open(outletFileName,'w')
outletFile.write('%i \n'%len(onl))
for i in range(len(onl)):
    nd = int(onl[i]); nd = nd;#make onl node numbers 0-based
    outletFile.write('%i \n'% nd) 
outletFile.close()


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
p_conv_fact = (((l_conv_fact/t_conv_fact)**2)*(1./3.))*(l_conv_fact**3)

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
    
    params.close()
    
else:
    print 'Run aborted.  Better luck next time!'


# ------------------------------------------------------------------------------------------
# add the partitioning work here; just use metis; if you can make it work with metis,
# it also should "just work" for geometric partitioning schemes

if lattice_type == 'D3Q15':
   lat = pp.D3Q15Lattice(int(Nx),int(Ny),int(Nz))
elif lattice_type == 'D3Q19':
   lat = pp.D3Q19Lattice(int(Nx),int(Ny),int(Nz))
else:
   lat = pp.D3Q27Lattice(int(Nx),int(Ny),int(Nz))


#lat15 = pp.D3Q15Lattice(int(Nx),int(Ny),int(Nz))
print "initializing the adjacency list"
lat.initialize_adjDict();
print "creating %s partition for %d processes" % (partition_style, numProcs)
lat.set_Partition(numParts= numProcs, style = partition_style)
lat.compute_cutSize()
print "cut size for %s partition = %g" % (partition_style, lat.get_cutSize())
#print "writing vtk file for %s partition" % partition_style
#partition_vtk_filename = "partition_%s.vtk" % partition_style
#lat.partition.write_vtk(partition_vtk_filename)
print "writing %s partition to disk" % partition_style
lat.partition.write_partition()
