import FluidChannel as fc
import numpy as np

# overall channel dimensions
aLx_p = 1.0
aLy_p = 0.75
aLz_p = 3
aNdivs = 5

# twin jet parameters
y1 = 0.21;
a = 0.05
S = 0.25
x1 = 0.25
W = 0.5
L = 2.

# create the obstruction object
TJ = fc.TwinJet(y1,x1,W,a,S,L)

# create the fluid channel object
myChan = fc.FluidChannel(Lx_p = aLx_p,Ly_p = aLy_p,Lz_p = aLz_p,
                         N_divs = aNdivs,obst = TJ);
                        
# write the mat file
myChan.write_mat_file('twinJet');

# write vtk of boundary conditions so you can visualize them
myChan.write_bc_vtk();