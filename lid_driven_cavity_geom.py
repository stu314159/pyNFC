# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 20:19:51 2017

@author: stu
"""

import FluidChannel as fc
import argparse

parser = argparse.ArgumentParser(prog='lid_driven_cavity_geom.py',
                                 description='create geometry files for lid-driven cavity problem')
                                 
parser.add_argument('nDivs',type=int)

# parse input arguments
args = parser.parse_args()
aN_divs = args.nDivs

aLx_p = 0.203 # meters
aLy_p = 0.240 # meters
aLz_p = 1.120 - 0.914 #meters

ld_cav = fc.LidDrivenCavity(Lx_p = aLx_p,Ly_p = aLy_p,
                            Lz_p = aLz_p,N_divs = aN_divs)
                            
ld_cav.write_bc_vtk()
ld_cav.write_mat_file('LDC_geom_test')
