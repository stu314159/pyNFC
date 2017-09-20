#!/usr/bin/env python
"""
geometry definition file for the gridChan geometry
"""
import sys
sys.path.insert(1,'.')

import FluidChannel as fc
import argparse

parser = argparse.ArgumentParser(prog='gridChan_geom.py',
                                 description='create geometry files gridChan problem')
                                 
parser.add_argument('nDivs',type=int)

# parse input arguments
args = parser.parse_args()
aN_divs = args.nDivs

# overall channel dimensions
aLx_p = 0.40 # meters
aLy_p = 0.40 # meters
aLz_p = 0.45 # meters (change to 2.15 m for actual LBM calcs)

# thickness of grids (m)
xT = 0.01 # x-thickness of vertical grids
yT = 0.01 # y-thickness of horizontal grids
zT = 0.01 # z-thickness of both horizontal and vertical grids

xPitch = 0.04 # x-pitch of vertical grids
yPitch = 0.04 # y-pitch of horizontal grids

gridZ = 0.35 # z-location of center of grid structure

hX = 0.20  # x-dimension of the circular hole in the grid
hY = 0.20  # y-dimension of the circular hole in the grid
hD = 0.043 # diameter of the circular hole in the grid

gridObst = fc.GridObst(gridZ,xT,yT,zT,xPitch,yPitch,hX,hY,hD);

gridChan = fc.FluidChannel(Lx_p = aLx_p, Ly_p = aLy_p, \
                           Lz_p = aLz_p, obst=gridObst, \
                           N_divs = aN_divs);

gridChan.write_mat_file('gridChan')

# visualize boundary condition data
gridChan.write_bc_vtk()

# set reference pressure near the end of the channel
gridChan.set_pRef_indx(aLx_p/2.,aLy_p/2.,0.98*aLz_p)