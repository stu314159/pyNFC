#!/p/home/sblair/anaconda2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 14:23:52 2017

@author: stu
"""

import sys
sys.path.insert(1,'.')

import pyPartition as pp
#from pymetis import part_graph #<-- requires that the PrgEnv-intel module be selected
import numpy as np
import scipy.io
import math
import argparse

parser = argparse.ArgumentParser(prog='pyNFC_partition.py',
                                 description='lattice partitioning script for pyNFC')

parser.add_argument('geom_filename',type=str)
parser.add_argument('lattice_type',type=str)
parser.add_argument('partition_style',type=str)
parser.add_argument('numProcs',type=int)

# parse input arguments
args = parser.parse_args()

# assign to required variables
geom_filename = args.geom_filename
lattice_type = args.lattice_type
partition_style = args.partition_style
numProcs = args.numProcs

geom_input = scipy.io.loadmat(geom_filename)
# overall domain dimensions
Lx_p = float(geom_input['Lx_p'])
Ly_p = float(geom_input['Ly_p'])
Lz_p = float(geom_input['Lz_p'])
Lo = float(geom_input['Lo'])
Ny_divs = int(geom_input['Ny_divs'])
rho_p = float(geom_input['rho_p'])
nu_p = float(geom_input['nu_p'])
#snl = list((geom_input['snl']).flatten())
#inl = list((geom_input['inl']).flatten()) # must be inlet on Z-min
#onl = list((geom_input['onl']).flatten()) # must be outlet on Z-max
ndType = list((geom_input['ndType']).flatten())


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

#print 'There are %d nodes in the solid node list'%len(snl)
#print 'Writing those nodes to file'
## now write this obstList to file.
#obstFilename = 'snl.lbm'
#obstFile = open(obstFilename,'w')
#obstFile.write('%i \n'%len(snl))
#for i in range(len(snl)):
#    nd = int(snl[i]); nd=nd;# make snl node numbers 0-based
#    obstFile.write('%i \n'% nd) 
#obstFile.close()
#
#print 'There are %d nodes in the inlet node list'%len(inl)
#print 'Writing those nodes to file'
#inletFileName = 'inl.lbm'
#inletFile = open(inletFileName,'w')
#inletFile.write('%i \n'%len(inl))
#for i in range(len(inl)):
#    nd = int(inl[i]); nd = nd;#make inl node numbers 0-based
#    inletFile.write('%i \n'% nd) 
#inletFile.close()
#
#print 'There are %d nodes in the outlet node list'%len(onl)
#print 'Writing those nodes to file'
#outletFileName = 'onl.lbm'
#outletFile = open(outletFileName,'w')
#outletFile.write('%i \n'%len(onl))
#for i in range(len(onl)):
#    nd = int(onl[i]); nd = nd;#make onl node numbers 0-based
#    outletFile.write('%i \n'% nd) 
#outletFile.close()

print 'There are %d nodes listed in ndType'%len(ndType)
print 'Writing those to file'
ndTypeFileName = 'ndType.lbm'
ndTypeFile = open(ndTypeFileName,'w')
for i in range(len(ndType)):
    nT = int(ndType[i]);
    ndTypeFile.write('%i \n'%nT)
ndTypeFile.close()

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
