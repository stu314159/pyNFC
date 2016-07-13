import FluidChannel as fc
import numpy as np

#ocStu = fc.FluidChannel(Lx_p = 2., Ly_p = 2., Lz_p = 8.,
#                        N_divs = 19,
#                        obst = fc.EllipticalScourPit(1.0,4.,0.1))
#ocStu.write_mat_file()

# geometric parameters for the spherical portion of the golf ball
d_golf_ball = 0.0427 # meters
aLx_p = 0.16 
aLy_p = 0.16
aLz_p = 0.3
aN_divs = 80

# construct the basic sphere
print 'Constructing the channel with smooth sphere'
sphereB = fc.SphereObstruction(d_golf_ball/2., aLx_p/2., aLy_p/2., aLz_p/2.)
sphereChannel = fc.FluidChannel(Lx_p = aLx_p,Ly_p = aLy_p, Lz_p = aLz_p,
                                N_divs = aN_divs, obst = sphereB)
sphereChannel.write_mat_file('sphere_n80')


# construct channel with the golf ball
print 'constructing the channel with the golf ball'


d_dimp = d_golf_ball/8.

# set the dimply 75% of its radius from surface
rd_dimp = (d_golf_ball/2.)+(d_dimp/2.)*0.65 
circ_golf_ball = np.pi*d_golf_ball;
N_dimp_max = np.floor(circ_golf_ball/d_dimp)

N_e = int(N_dimp_max*0.6)
N_a = 2*N_e-1 # for now, assume they are the same density in each direction
# since the azimuths will go all the way around and elevations
# only half way; use this formulation

# construct the golf ball
golfB = fc.GolfBall(sphereB,d_dimp,rd_dimp,N_e,N_a)

golfChannel = fc.FluidChannel(Lx_p = aLx_p,Ly_p = aLy_p,Lz_p = aLz_p,
                               N_divs = aN_divs, obst = golfB)

#golfChannel.write_bc_vtk()
golfChannel.write_mat_file('golfBall_n80')

