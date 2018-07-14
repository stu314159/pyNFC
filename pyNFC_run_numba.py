#!/usr/bin/env python
##!/home/users/sblair/anaconda2/bin/python
"""
 testing implementation of pyNFC
"""

from mpi4py import MPI
import time

#from vtkHelper import saveVelocityAndPressureVTK_binary as writeVTK
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
import sys
sys.path.insert(1,'.')
import pyNFC_numba
# Information about the LBM run that produced the data - I should get this from params.lbm
# Read data from params.lbm
input_file_name = 'params.lbm'
input_data = open(input_file_name,'r')
lattice_type = str(input_data.readline().rstrip()) # perhaps I can actually use this as intended now
dynamics = int(input_data.readline())
Num_ts = int(input_data.readline())
ts_rep_freq = int(input_data.readline())
Warmup_ts = int(input_data.readline())
plot_freq = int(input_data.readline())
Cs = float(input_data.readline())
rho_lbm = float(input_data.readline())
u_lbm = float(input_data.readline())
omega = float(input_data.readline())
Nx = int(input_data.readline())
Ny = int(input_data.readline())
Nz = int(input_data.readline())
Restart_flag = int(input_data.readline())
TimeAvg_flag = int(input_data.readline())
Lx_p = float(input_data.readline())
Ly_p = float(input_data.readline())
Lz_p = float(input_data.readline())
t_conv_fact = float(input_data.readline())
l_conv_fact = float(input_data.readline())
p_conv_fact = float(input_data.readline())
pRef_ind = int(input_data.readline())
SubspaceData_flag = int(input_data.readline())
input_data.close()

#print "process %d of %d says hello!" % (rank, size)

# each process initialize their partition:

myPart = pyNFC_numba.NFC_LBM_partition(rank,size,comm,Nx,Ny,Nz,rho_lbm,u_lbm,dynamics,omega,Cs,lattice_type)

myPart.write_node_sorting() # should this be done in the constructor?

# load restart data if applicable
if Restart_flag == 1:
    myPart.load_restart_data()
    
if TimeAvg_flag == 1:
    myPart.initialize_timeAvg()
    
if SubspaceData_flag == 1:
    myPart.allocate_subspace_data_arrays(Num_ts)
    

# if there are subset nodes and if subset data is to be stored,
# initialize the data structure that will store the subset data 
# and write it to disk in an appropriate data structure.
    
# NFC_LBM_partition should own this data structure.
    
# do some time stepping
#numTs = 10
#plot_freq = 5;
if rank == 0:
    if ((Restart_flag==1) and (rank==0)):
        print "Loading restart data."
    time1 = time.time()

for ts in range(Num_ts):
    
    if (((ts+1)%ts_rep_freq)==0):
    	if rank==0:
            print "executing time step %d" % (ts+1)
    isEven = (ts%2 == 0)
    myPart.take_LBM_timestep(isEven)
    
    if (SubspaceData_flag == 1):
        myPart.record_subspace_data(ts)

    if ((ts % plot_freq == 0) and (ts > Warmup_ts)):
        myPart.write_data(isEven)
    

if rank == 0:
    time2 = time.time()
    ex_time = time2 - time1
    print "approximate execution time is %g seconds" % (ex_time)
    numLP = Nx*Ny*Nz
    LPUs = numLP*Num_ts/(ex_time)
    print "approximate LPU/sec = %g " % LPUs
    
if SubspaceData_flag == 1:
    if rank == 0:
        print "Writing subspace data arrays to binary data files"
        
    myPart.write_subspace_data();
    myPart.write_ss_node_sorting();
    
if TimeAvg_flag == 1:
    myPart.write_timeAvg();
