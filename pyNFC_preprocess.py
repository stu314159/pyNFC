#pyNFC_preprocess.py
"""
   pre-processing steps for pyNFC code.  Basic geometry definition
   and partitioning information created and input file created prior to
   parallel simulation run

"""

import FluidChannel as fc
import pyNFC 
import pyPartition as pp



# define geometry - channel size, obstruction and FluidChannel properties
Lx_p = 1.; Ly_p = 1.; Lz_p = 8.; N_divs = 8;
print "getting a spherical obstruction"
obstacle = fc.SphereObstruction(r=Lx_p/8., x_c = Lx_p/2.,
                                y_c = Ly_p/2., z_c = Lz_p/2.)
print "creating the fluid channel"
chanGeom = fc.FluidChannel(Lx_p = Lx_p, Ly_p = Ly_p, Lz_p = Lz_p,
                           N_divs = N_divs, obst = obstacle) 
# set desired number of partitions:
numPartitions = 16; numTrials = 2000;

print "creating partitions"
# call parititioning libraries to obtain a decent geometric partition
myPart = pp.Partitioner(Nx = int(chanGeom.Nx), Ny = int(chanGeom.Ny), Nz = int(chanGeom.Nz),
                        numParts = numPartitions, numTrials = numTrials)
print "writing partitions to file"
myPart.write_vtk()

