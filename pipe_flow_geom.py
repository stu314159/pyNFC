# -*- coding: utf-8 -*-
"""

Pipe flow problem
Created on Thu Aug  3 11:32:57 2017

@author: stu
"""
import FluidChannel as fc
import argparse

parser = argparse.ArgumentParser(prog='pipe_flow_geom.py',
                                 description='create geometry files for pipe flow problem')
                                 
parser.add_argument('nDivs',type=int)

# parse input arguments
args = parser.parse_args()

#overall channel dimensions
aLx_p = 1.1
aLy_p = 1.1
aLz_p = 10.0
aNdivs = args.nDivs

x_c = aLx_p/2.;
y_c = aLy_p/2.;
pDiam = min([aLx_p,aLy_p]) - 0.1

myObst = fc.StraightPipe(x_c,y_c,pDiam);

myChan = fc.FluidChannel(Lx_p=aLx_p,Ly_p=aLy_p,Lz_p=aLz_p,obst=myObst,
                         N_divs=aNdivs)                         

# write the mat file
myChan.write_mat_file('pipe_flow');

# write vtk of boundary conditions so you can visualize them
#myChan.write_bc_vtk();

myChan.set_pRef_indx(aLx_p/2.,aLy_p/2.,0.95*aLz_p)
