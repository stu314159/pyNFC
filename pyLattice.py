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

         input:
           f - local particle density distribution

         outputs:
           rho - macroscopic density
           ux, uy, uz - macroscopic velocity components.
        """
        rho = np.sum(f)
        ux = np.dot(f,self.ex)/rho;
        uy = np.dot(f,self.ey)/rho;
        uz = np.dot(f,self.ez)/rho;
        return rho, ux, uy, uz  

    def compute_equilibrium(self,f,rho,u):
        """
           given f, rho and u (vector - all macro velocities)
           return f_eq

        """
        numSpd = self.get_numSpd();
        ux = u[0]; uy = u[1]; uz = u[2];
        f_eq = np.zeros(numSpd,dtype = np.float32)
        for spd in numSpd:
             cu = 3.*(self.ex[spd]*ux + self.ey[spd]*uy + self.ez[spd]*uz)
             f_eq[spd] = self.w[spd]*rho*(1. + cu + 0.5*(cu*cu) - 
                           3./2.*(ux*ux + uy*uy + uz*uz)) 

        return f_eq[:]   


    def compute_fOut(self,fIn,ndType,omega,Cs=0.,u_bc=0.,rho_bc=1.):
        """
           input:
            fIn - incoming particle density distribution
            ndType - [0 | 1 | 2 | 3] 
                     0 = interior fluid node
                     1 = solid node
                     2 = inlet velocity node
                     3 = outlet pressure node
            omega - relaxation parameter
            Cs - (optional) Smagorinsky turbulence parameter
            u_bc - (optional) inlet velocity 
            rho_bc - (optional) outlet pressure
            
          output:
           fOut - outgoing particle density distribution for streaming

        """
        fluidNode = False; solidNode = False; inletNode = False; outletNode = False;
        if ndType == 0:
            fluidNode = True;
        elif ndType == 1:
            solidNode = True;
        elif ndType == 2:
            inletNode = True;
        elif ndType == 3:
            outletNode = True;
        else:
            raise ValueError("ndType must be one of 0, 1, 2, or 3")

        rho, ux, uy, uz = self.compute_macroscopic_data(fIn)
        if ndType == 2:
            rho = self.set_inlet_velocity_bc_macro(fIn,uz)
        elif ndType == 3:
            uz = self.set_outlet_pressure_bc_macro(fIn,rho)

        if ndType != 1: #solid nodes do not need fEq
            fEq = self.compute_equilibrium(fIn,rho,[ux,uy,uz])

        
        if ((ndType == 2) or (ndType == 3)):
            pass # implement soon



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
        rho = (1./(1. - uz))*(2.*(f[6]+f[11]+f[12]+f[13]+f[14])+(f[0]+f[1]+f[2]+f[3]+f[4]))

        return rho
 

    def set_outlet_density_bc_macro(self,f,rho): 
        """
          compute macroscopic uz for density outlet
          bc using Regularized BC methods
        """
        uz = -1. + (1./rho)*(2.*(f[6]+f[11]+f[12]+f[13]+f[14])+(f[0]+f[1]+f[2]+f[3]+f[4]))
    
        return uz

    def bounceBack_inletBoundary_micro(self,f,fEq):
        """
          input:
             f and fEq

          output:
             f with micro velocities for speeds into the
             domain (unknown) adjusted by "bouncing back"
             the non-equilibrium component of the known
             speeds in opposite direction:


          for D3Q15, unknown speeds on (low-z) inlet: 5, 7, 8, 9, 10
                 corresponding bounce-back speeds: 6, 14, 13, 12, 11

        """
        sp = [5,7,8,9,10]; bbSp = [6,14,13,12,11];
        f[sp] += f[bbSp] - fEq[bbSp];


    def bounceBack_outletBoundary_micro(self,f,fEq):
        pass #implement soon

    def regularize_boundary_nodes(self,f,fEq):
        pass #implement soon




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
        return rho

    def set_outlet_density_bc_macro(self,f,rho): 
        """
          compute macroscopic uz for density outlet
          bc using Regularized BC methods
        """
        uz = -1. + (1./rho)*(2.*(f[6]+f[13]+f[14]+f[17]+f[18])+(f[0]+f[1]+f[2]+f[3]+f[4]+f[7]+f[8]+f[9]+f[10]))

        return uz

    def bounceBack_inletBoundary_micro(self,f,fEq):
        """
          input:
             f and fEq

          output:
             f with micro velocities for speeds into the
             domain (unknown) adjusted by "bouncing back"
             the non-equilibrium component of the known
             speeds in opposite direction:


          for D3Q19, unknown speeds on (low-z) inlet: 5,11,12,15,16
                 corresponding bounce-back speeds: 6,14,13,18,17

        """
        sp = [5,11,12,15,16]; bbSp = [6,14,13,18,17];
        f[sp] += f[bbSp] - fEq[bbSp];

    def bounceBack_outletBoundary_micro(self,f,fEq):
        pass #implement soon

    def regularize_boundary_nodes(self,f,fEq):
        pass #implement soon



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
        rho = (1./(1. - uz))*(2.*(f[3]+f[6]+f[8]+f[10]+f[12]+f[20]+f[22]+f[24]+f[26])+ 
                             (f[0]+f[1]+f[2]+f[4]+f[5]+f[14]+f[15]+f[17]+f[18]))
        return rho


    def set_outlet_density_bc_macro(self,f,rho): 
        """
          compute macroscopic uz for density outlet
          bc using Regularized BC methods
        """
        uz = -1. + (1./rho)*(2.*(f[3]+f[6]+f[8]+f[10]+f[12]+f[20]+f[22]+f[24]+f[26])+ 
                             (f[0]+f[1]+f[2]+f[4]+f[5]+f[14]+f[15]+f[17]+f[18]))
        return uz


    def bounceBack_inletBoundary_micro(self,f,fEq):
         """
          input:
             f and fEq

          output:
             f with micro velocities for speeds into the
             domain (unknown) adjusted by "bouncing back"
             the non-equilibrium component of the known
             speeds in opposite direction:


          for D3Q27, unknown speeds on (low-z) inlet: 7,9,11,13,16,19,21,23,25
                          bbSpd = 20,22,24,26,3,6,8,10,12
                 

        """
        sp = [7,9,11,13,16,19,21,23,25]; bbSp = [20,22,24,26,3,6,8,10,12];
        f[sp] += f[bbSp] - fEq[bbSp];


    def bounceBack_outletBoundary_micro(self,f,fEq):
        pass #implement soon

    def regularize_boundary_nodes(self,f,fEq):
        pass #implement soon



if __name__=="__main__":
    """  
      put testing code here
    """

from sympy import *

f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f13 = symbols("f0 f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14")
print "f0 + f1 = " + str(f0 + f1)
