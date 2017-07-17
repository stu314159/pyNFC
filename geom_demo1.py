import FluidChannel as fc
import numpy as np

#overall channel dimensions
aLx_p = 1.0
aLy_p = 1.0
aLz_p = 5.0
aNdivs = 21

#sphere position
a_x = 0.5
a_y = 0.5
a_z = 2.0
r = 0.2

#create the obstruction object
myObst = fc.SphereObstruction(r,a_x,a_y,a_z);

#creat the fluid channel object
myChan = fc.FluidChannel(Lx_p = aLx_p,Ly_p = aLy_p,Lz_p = aLz_p,
                         N_divs = aNdivs,obst = myObst);

# write the mat file
myChan.write_mat_file('demo1');

# write vtk of boundary conditions so you can visualize them
myChan.write_bc_vtk();

