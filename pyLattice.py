# pyLattice.py
"""
class module for the Lattice class to be used with pyNFC


"""
import numpy as np

class Lattice(object):
    """
       define the layout and adjacency of the LBM lattice
    """
    def __init__(self,Nx,Ny,Nz):
        """
            basic constructor
            Nx - number of lattice points in the x-direction
            Ny - number of lattice points in the y-direction
            Nz - number of lattice points in the z-direction
            
        """
        self.Nx = Nx; self.Ny = Ny; self.Nz = Nz
        self.ex = []; self.ey = []; self.ez = []
        self.bbSpd = []; self.w = []; 
        

    def get_dims(self):
        return [self.Nx, self.Ny, self.Nz]

    def get_ex(self):
        return self.ex[:]


    def get_ey(self):
        return self.ey[:]

    def get_ez(self):
        return self.ez[:]

    def get_numSpd(self):
        return len(self.ex)

    def get_bbSpd(self):
        return self.bbSpd[:]

    def get_w(self):
        return self.w[:]

    def compute_macroscopic_data(self,f):
        """
          given microscopic density distribution values,
          compute and return macroscopic density and velocity
        """
        rho = np.sum(f)
        ux = np.dot(f,self.ex)/rho;
        uy = np.dot(f,self.ey)/rho;
        uz = np.dot(f,self.ez)/rho;
        return rho, ux, uy, uz  
    

class D3Q15Lattice(Lattice):
    """
      D3Q15 Lattice
    """
    def __init__(self,Nx,Ny,Nz):
        """
          D3Q15 Lattice
        """
        super(D3Q15Lattice,self).__init__(Nx,Ny,Nz)
        self.ex =  [0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1]
        self.ey = [0,0,0,1,-1,0,0,1,1,-1,-1,1,1,-1,-1]
        self.ez = [0,0,0,0,0,1,-1,1,1,1,1,-1,-1,-1,-1]

        self.bbSpd = [0,2,1,4,3,6,5,14,13,12,11,10,9,8,7]
        self.w = [2./9.,1./9.,1./9,1./9.,1./9.,1./9.,1./9.,
	       1./72.,1./72.,1./72.,1./72.,
	       1./72.,1./72.,1./72.,1./72.]

    def set_inlet_velocity_bc_macro(self,f,uz): # not too flexible, but it is what NFC does (one thing at a time)
        """
          compute macroscopic density for velocity inlet
          bc using Regularized BC methods
        """

        #rho = (1./(1.-uz))*(2.0*(f6+f11+f12+f13+f14)+(f0+f1+f2+f3+f4)); <-- C code
        rho = (1./(1. - uz))*(2.*(f[6]+f[11]+f[12]+f[13]+f[14])+(f[0]+f[1]+f[2]+f[3]+f[4])
    


class D3Q19Lattice(Lattice):
    """
    """
    def __init__(self,Nx,Ny,Nz):
        super(D3Q19Lattice,self).__init__(Nx,Ny,Nz)
        self.ex =  [0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0]
        self.ey = [0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1]
        self.ez = [0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1]
        self.bbSpd = [0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15]
        self.w = [2./9.,1./9.,1./9,1./9.,1./9.,1./9.,1./9.,
	       1./72.,1./72.,1./72.,1./72.,
	       1./72.,1./72.,1./72.,1./72.]

    def set_inlet_velocity_bc_macro(self,f,uz):
        """
          compute macroscopic density for velocity inlet bc using
          Regularized BC methods
        """
        #rho = (1./(1.-uz))*(2.0*(f6+f13+f14+f17+f18)+(f0+f1+f2+f3+f4+f7+f8+f9+f10));
        rho = (1./(1.-uz))*(2.*(f[6]+f[13]+f[14]+f[17]+f[18])+(f[0]+f[1]+f[2]+f[3]+f[4]+f[7]+f[8]+f[9]+f[10]))

   
class D3Q27Lattice(Lattice):
    """
    """
    def __init__(self,Nx,Ny,Nz):
        super(D3Q27Lattice,self).__init__(Nx,Ny,Nz)
        self.ex = [0,-1,0,0,-1,-1,-1,-1,0,0,-1,-1,-1,-1,1,0,0,1,1,1,1,0,0,1,1,1,1]
        self.ey = [0,0,-1,0,-1,1,0,0,-1,-1,-1,-1,1,1,0,1,0,1,-1,0,0,1,1,1,1,-1,-1]
        self.ez = [0,0,0,-1,0,0,-1,1,-1,1,-1,1,-1,1,0,0,1,0,0,1,-1,1,-1,1,-1,1,-1]
        self.bbSpd = [0,14,15,16,17,18,19,20,21,22,23,24,25,26,
	      1,2,3,4,5,6,7,8,9,10,11,12,13]
        self.w = [8./27.,2./27.,2./27.,2./27.,1./54.,1./54.,1./54.,1./54.,1./54.,
	       1./54.,1./216.,1./216,1./216.,1./216.,2./27.,2./27.,
	       2./27.,1./54.,1./54.,1./54.,1./54.,1./54.,
		1./54.,1./216.,1./216,1./216.,1./216.]

    def set_inlet_velocity_bc_macro(self,f,uz):
        pass  # update this when ready
