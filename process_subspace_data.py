#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 08:51:44 2018

@author: stu
"""

import sys
sys.path.insert(1,'.')

#import math
#import os
import numpy as np
import scipy.io
from hdf5Helper import writeH5, writeXdmf

# need to read subspace partition data and determine the number of time steps
input_file_name = 'params.lbm'
input_data = open(input_file_name,'r')
lattice_type = str(input_data.readline().rstrip()) #<- not needed
dynamics = int(input_data.readline()) #<- not needed
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
pRef_idx = int(input_data.readline())
input_data.close()

u_conv_fact = t_conv_fact/l_conv_fact;
nnodes = Nx*Ny*Nz
x = np.linspace(0.,Lx_p,Nx).astype(np.float64);
y = np.linspace(0.,Ly_p,Ny).astype(np.float64);
z = np.linspace(0.,Lz_p,Nz).astype(np.float64);
numEl = Nx*Ny*Nz
Y,Z,X = np.meshgrid(y,z,x);
XX = np.reshape(X,numEl)
YY = np.reshape(Y,numEl)
ZZ = np.reshape(Z,numEl)
dx = x[1] - x[0]

# need to get more information on the structure of the subset data file
ss_filename = 'part_ss_sizes.mat' 
ss_data = scipy.io.loadmat(ss_filename)
part_ss_sizes = ss_data['part_ss_sizes'] # lists the number of ss nodes in each partition
#num_ss_partitions = len(part_ss_sizes); # number of partitions used in the simulation
total_ss_nodes = sum(part_ss_sizes.flatten());

part_ss_sizes = np.array(part_ss_sizes,dtype=np.int32).reshape(4,);
num_ss_partitions = len(part_ss_sizes)
ss_offsets = np.cumsum(part_ss_sizes)



# I want to put these into 4 numpy arrays.  Each array will be of size
# equal to the total number of subspace nodes x num_ts
ss_ux_raw = np.fromfile('ss_ux.b_dat',dtype=np.float32);
#ss_ux_i = np.reshape(ss_ux_raw, (Num_ts,int(total_ss_nodes)))

ss_uy_raw = np.fromfile('ss_uy.b_dat',dtype=np.float32);
#ss_uy_i = np.reshape(ss_uy_raw,(Num_ts,int(total_ss_nodes)));

ss_uz_raw = np.fromfile('ss_uz.b_dat',dtype=np.float32);
#ss_uz_i = np.reshape(ss_uz_raw,(Num_ts,int(total_ss_nodes)));

ss_rho_raw = np.fromfile('ss_rho.b_dat',dtype=np.float32);
#ss_rho_i = np.reshape(ss_rho_raw,(Num_ts,int(total_ss_nodes)));

ss_ux_i = np.empty((Num_ts,int(total_ss_nodes)),dtype=np.float32);
ss_uy_i = np.empty_like(ss_ux_i)
ss_uz_i = np.empty_like(ss_ux_i)
ss_rho_i = np.empty_like(ss_ux_i)

for i in range(num_ss_partitions):
    if i==0:
        l_lp = 0
    else:
        l_lp = ss_offsets[i-1]
    
    u_lp = ss_offsets[i]
    
    l_idx = l_lp*Num_ts;
    u_idx = u_lp*Num_ts;
    
    ss_ux_i[:,l_lp:u_lp]=ss_ux_raw[l_idx:u_idx].reshape(Num_ts,int(part_ss_sizes[i]))
    ss_uy_i[:,l_lp:u_lp]=ss_uy_raw[l_idx:u_idx].reshape(Num_ts,int(part_ss_sizes[i]))
    ss_uz_i[:,l_lp:u_lp]=ss_uz_raw[l_idx:u_idx].reshape(Num_ts,int(part_ss_sizes[i]))


# convert units
ss_ux_i /= u_conv_fact;
ss_uy_i /= u_conv_fact;
ss_uz_i /= u_conv_fact;
ss_rho_i *= p_conv_fact;

# global node numbers of ss data in the order presented
g_order_map = np.fromfile('ss_ordering.b_dat',dtype=np.int32).astype(np.int32)

# re-order data so that node numbers are always increasing.
ordered_idx = np.argsort(g_order_map);

ss_ux = ss_ux_i[:,ordered_idx];
ss_uy = ss_uy_i[:,ordered_idx];
ss_uz = ss_uz_i[:,ordered_idx];
ss_rho = ss_rho_i[:,ordered_idx];

# each row of the above arrays consitutes the data for the subspace on a single time step.
# need Ny_ss and Nz_ss so we know the structure of this data. 
x_ss = XX[g_order_map[ordered_idx]];
y_ss = YY[g_order_map[ordered_idx]];
z_ss = ZZ[g_order_map[ordered_idx]];

X_ss = np.min(x_ss); # is constant

Y_ss_min = np.min(y_ss);
Y_ss_max = np.max(y_ss);
Ny_ss = (Y_ss_max - Y_ss_min)/dx + 1;

Z_ss_min = np.min(z_ss);
Z_ss_max = np.max(z_ss);
Nz_ss = (Z_ss_max - Z_ss_min)/dx + 1;

# re-shape the arrays into 3-D arrays: num_ts x Nz x Ny (?)
Nz_ss = int(Nz_ss)
Ny_ss = int(Ny_ss)

ss_ux_3 = np.reshape(ss_ux,(Num_ts,Nz_ss,Ny_ss))
ss_uy_3 = np.reshape(ss_uy,(Num_ts,Nz_ss,Ny_ss))
ss_uz_3 = np.reshape(ss_uz,(Num_ts,Nz_ss,Ny_ss))
ss_rho_3 = np.reshape(ss_rho,(Num_ts,Nz_ss,Ny_ss))

# re-order the data so it is amenable with Matlab/Octave
ss_ux_3 = np.array(ss_ux_3,order='F')
ss_uy_3 = np.array(ss_uy_3,order='F')
ss_uz_3 = np.array(ss_uz_3,order='F')
ss_rho_3 = np.array(ss_rho_3,order='F')

# save to a *.mat file
scipy.io.savemat('ss_data.mat',{'ss_ux':ss_ux_3, 'ss_uy':ss_uy_3, 
                                'ss_uz':ss_uz_3, 'ss_rho':ss_rho_3,
                                'dx':dx, 'dt':(1./t_conv_fact)})

velmag = np.sqrt(ss_ux**2 + ss_uy**2 + ss_uz**2);

# conver the time steps into a series of *.h5 files and associated *.xmf
dims = (Nz_ss,Ny_ss,1)
for ts in range(Num_ts):
    print 'processing subspace data set %d'%ts
    h5_file = 'ss_out'+ str(ts) + '.h5'
    xmf_file = 'ss_data' + str(ts) + '.xmf'
    
    ss_ux1 = ss_ux[ts,:]
    ss_uy1 = ss_uy[ts,:]
    ss_uz1 = ss_uz[ts,:]
    ss_pressure = ss_rho[ts,:]
    ss_pressure = ss_pressure.reshape(dims)
    velmag_1 = velmag[ts,:]
    writeH5(ss_pressure,ss_ux1,ss_uy1,ss_uz1,velmag_1,h5_file)
    writeXdmf(dims,dx,xmf_file,h5_file)
