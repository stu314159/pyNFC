#!/usr/bin/env python
"""
 testing implementation of pyNFC
"""

from mpi4py import MPI
import numpy as np
import math


#from vtkHelper import saveVelocityAndPressureVTK_binary as writeVTK
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
import sys
sys.path.insert(1,'.')
from vtkHelper import saveVelocityAndPressureVTK_binary as writeVTK
import pyNFC
import FluidChannel as fc
# Information about the LBM run that produced the data - I should get this from params.lbm
# Read data from params.lbm
input_file_name = 'params.lbm'
input_data = open(input_file_name,'r')
latticeType = int(input_data.readline())
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
Lx_p = float(input_data.readline())
Ly_p = float(input_data.readline())
Lz_p = float(input_data.readline())
t_conv_fact = float(input_data.readline())
l_conv_fact = float(input_data.readline())
p_conv_fact = float(input_data.readline())
input_data.close()

# each process initialize their partition:
lattice_type = 'D3Q27'
myPart = pyNFC.NFC_LBM_partition(rank,size,Nx,Ny,Nz,lattice_type)
