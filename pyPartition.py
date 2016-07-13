#pyPartition.py
"""
provide class libraries for partitioning
goal of this code is to provide a reasonable geometric partition
to the pre-processing libraries.  The output will be a list of length Nx*Ny*Nz
containing the integer of which partition each lattice point lives.

"""

import numpy as np
import partition_suggestion as ps
import partition_compare as pc
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTK

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
        self.adjDict = {};

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

    def get_nnodes(self):
        return self.Nx*self.Ny*self.Nz

    def initialize_adjDict(self):  # just playing around here
        if len(self.adjDict.keys()) != 0:
            self.adjDict = pc.set_adjacency(self.Nx,self.Ny,self.Nz,self.ex,self.ey,self.ez)

    def add_Partition(self):
        self.partition = Partitioner(self.Nx, self.Ny, self.Nz



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

   
class D3Q27Lattice(Lattice):
    """
    """
    def __init__(self,Nx,Ny,Nz):
        super(D3Q27Lattice,self).__init__(Nx,Ny,Nz)
        self.ex =  [0,-1,0,0,-1,-1,-1,-1,0,0,-1,-1,-1,-1,1,0,0,1,1,1,1,0,0,1,1,1,1]
        self.ey = [0,0,-1,0,-1,1,0,0,-1,-1,-1,-1,1,1,0,1,0,1,-1,0,0,1,1,1,1,-1,-1]
        self.ez = [0,0,0,-1,0,0,-1,1,-1,1,-1,1,-1,1,0,0,1,0,0,1,-1,1,-1,1,-1,1,-1]
        self.bbSpd = [0,14,15,16,17,18,19,20,21,22,23,24,25,26,
	      1,2,3,4,5,6,7,8,9,10,11,12,13]
        self.w = [8./27.,2./27.,2./27.,2./27.,1./54.,1./54.,1./54.,1./54.,1./54.,
	       1./54.,1./216.,1./216,1./216.,1./216.,2./27.,2./27.,
	       2./27.,1./54.,1./54.,1./54.,1./54.,1./54.,
		1./54.,1./216.,1./216,1./216.,1./216.]

class Partitioner:
    """
     the class that will do the work to select and obtain a partition
    """

    def __init__(self,Nx,Ny,Nz,numParts,numTrials = 2000):
        """
          Nx - number of lattice points in the x-direction (int)
          Ny - number of lattice points in the y-direction (int)
          Nz - number of lattice points in the z-direction (int)
          numParts - number of partitions to form (int)
          numTrials - number of attempts that the randomized
                      partition advisor should use to find 
                      a good partitioning
                      
        """
        self.Nx = Nx; self.Ny = Ny; self.Nz = Nz
        self.numParts = numParts
        self.numTrials = numTrials
        print "calling part advisor"
        [self.px,self.py,self.pz] = ps.part_advisor(self.Nx,self.Ny,self.Nz,self.numParts,self.numTrials)
        print "setting geometric partition"
        self.part_vert = pc.set_geometric_partition(self.Nx, self.Ny, self.Nz,
                                                self.px, self.py, self.pz)



    def get_partition(self):
        """
          give access to partition
        """

        return self.part_vert[:]


    def get_partition_sizes(self):
        """
          give access to partition sizes
        """
        return [self.px, self.py, self.pz]

    def write_vtk(self):
        """
          write out a vtk file to allow visualization of the partitioning
        """
        dims = [self.Nx, self.Ny, self.Nz]
        origin = [0., 0., 0.]
        spacing = [0.1, 0.1, 0.1] #<-- no need for this to correspond to actual physical spacing
        writeVTK(self.part_vert,'partitions','partition_pa.vtk',dims,origin,spacing)



   
        


if __name__=="__main__":
    """
      put testing code here
    """
    lat = D3Q15Lattice(5,5,5)
    print "ex = ", lat.get_ex()
    print "num_spd = ", lat.get_numSpd()

    lat2 = D3Q19Lattice(5,5,5)
    print "num_spd D3Q19 = ", lat2.get_ex()
    print "num_spd D3Q19 = ", lat2.get_numSpd()
