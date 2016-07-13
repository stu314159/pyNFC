#pyPartition.py
"""
provide class libraries for partitioning
goal of this code is to provide a reasonable geometric partition
to the pre-processing libraries.  The output will be a list of length Nx*Ny*Nz
containing the integer of which partition each lattice point lives.

"""

import partition_suggestion as ps
import partition_compare as pc
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTK

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



   
        


if __name__=='main':
    """
      put testing code here
    """
