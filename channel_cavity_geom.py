#!/p/home/sblair/anaconda2/bin/python

import sys
sys.path.insert(1,'.')

import FluidChannel as fc
import argparse

parser = argparse.ArgumentParser(prog='channel_cavity_geom.py',
                                 description='create geometry files for cavity channel problem')
                                 
parser.add_argument('nDivs',type=int)

# parse input arguments
args = parser.parse_args()
aN_divs = args.nDivs

# overall channel dimensions
aLx_p = 0.203 # meters
aLy_p = 0.253 # meters
aLz_p = 2.00 # meters

# cavity parameters
cDepth = 0.240 # meters
cStart = 0.914 # meters
cEnd = 1.120 # meters

cavity = fc.ChannelCavity(cDepth,cStart,cEnd)
testChan = fc.FluidChannel(Lx_p = aLx_p, Ly_p = aLy_p, Lz_p = aLz_p,
                           N_divs = aN_divs, obst = cavity)
testChan.write_mat_file('ChanCavityTest')
#testChan.write_bc_vtk()

testChan.set_pRef_indx(aLx_p/2.,0.2465,0.98*aLz_p)

