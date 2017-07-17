#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 08:54:21 2017

@author: sblair
"""

import FluidChannel as fc

#overall channel dimensions
aLx_p = 6.4
aLy_p = 3.0
aLz_p = 14.0
aNdivs = 10

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
myChan.write_mat_file('wall_mounted_brick');

# write vtk of boundary conditions so you can visualize them
myChan.write_bc_vtk();