#pyNFC.py
"""
Implementation file for pyNFC library module -- under development

"""

import numpy as np
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTKpt
from mpi4py import MPI

class NFC_Halo_Data_Organizer(object):
    """
     collect and organize how halo data will be organized prior to communication
    """
    def __init__(self,ngb_rank):
        """
          constructor
        """
        self.ngb_rank = ngb_rank


class NFC_Part_Communicator(object):
    """
     class designed to handle the communication tasks for an NFC_LBM_partition
    """
    def __init__(self,rank,size,comm,comm_list_out):
        """
          rank - which MPI process
          size - number of MPI processes
          comm - MPI communicator
          comm_list - list of tuples containing all data to be exchanged
          

          Each partition needs to know:
            a) the list of MPI ranks it needs to exchange halo data with;
            b) for each exchange pair the # of data elements to be exchanged;
            
            Of course each process also has to manage the local/global node numbers
            for which each data element is bound as well as the associated speed.

            (for both incoming and outgoing data) as well as what order it will arrive in.
            by convention: 

                a) each partition will send data to neighboring partitions in order of 
                   increasing global lattice point number; and
                b) for lattice points receiving multiple speeds from the SAME neighbor partition,
                   the data will be provided by increasing spd number.

            how this is accomplished, and which classes accomplish this, is the design question
            yet to be answered.
          
        """
        self.rank = rank; self.size = size; self.comm = comm;
        self.comm_list = comm_list;
        self.lattice = lattice;
        




class NFC_LBM_partition(object):
    """
    each partition has:
         
    """
    def __init__(self,rank,size,comm,Nx,Ny,Nz,lattice_type='D3Q15'):
        """
          rank - MPI rank for this partition
          size - MPI size for MPI COMM WORLD
          lattice_type - ['D3Q15' | 'D3Q19' | 'D3Q27']
        """
        self.rank = rank; self.size = size; self.comm = comm # MPI variabes
        self.Nx = Nx; self.Ny = Ny; self.Nz = Nz; # all partitions need to know the global domain structure
        
        if lattice_type == 'D3Q15':
            self.lattice = D3Q15Lattice(self.Nx, self.Ny, self.Nz)
        elif lattice_type == 'D3Q19':
            self.lattice = D3Q19Lattice(self.Nx, self.Ny, self.Nz)
        else:
            self.lattice = D3Q27Lattice(self.Nx, self.Ny, self.Nz)

        self.numSpd = self.lattice.get_numSpd()
        #print "process %d of %d constructed %s lattice " % (rank,size,lattice_type)
        self.ex = np.array(self.lattice.get_ex(),dtype=np.int32);
        self.ey = np.array(self.lattice.get_ey(),dtype=np.int32);
        self.ez = np.array(self.lattice.get_ez(),dtype=np.int32);
        self.bb_Spd = np.array(self.lattice.get_bbSpd(),dtype=np.int32);
        self.load_parts()
        self.gen_adjacency() # adjacency list using global node numbers
        self.get_halo_nodes() # halo nodes are all global node numbers
        self.write_partition_vtk() # visualize each partition interior, boundary and halo


    def get_XYZ_index(self,g_nd): # this will depend upon a global geometry structure like a brick
        z = g_nd/(self.Nx*self.Ny)
        y = (g_nd - z*self.Nx*self.Ny)/self.Nx
        x = (g_nd - z*self.Nx*self.Ny - y*self.Nx)
        return (x,y,z)

    def get_gInd_XYZ(self,x,y,z): # this will give global index given x, y, z index
        return x+y*self.Nx + z*self.Nx*self.Ny

    def gen_adjacency(self):
        """
          create an adjacency matrix holding the global node numbers of 
          neighbors to all lattice points
        """
        self.adjacency = np.empty((self.num_local_nodes,self.numSpd),dtype=np.int32) # create the array
        
        # do this the bone-headed way:
        for k in self.global_to_local.keys():
            (x,y,z) = self.get_XYZ_index(k)
            for spd in range(self.numSpd):
                x_t = x + self.ex[spd]; x_t = x_t % self.Nx;
                y_t = y + self.ey[spd]; y_t = y_t % self.Ny;
                z_t = z + self.ez[spd]; z_t = z_t % self.Nz;
                self.adjacency[self.global_to_local[k],spd] = self.get_gInd_XYZ(x_t,y_t,z_t)
                
                


    def get_halo_nodes(self):
        """
         Identify global node numbers that must be incorporated onto the halo.
         Identify the stream-out speeds associated with each halo node and the partition
         each halo node belongs in.
         Also compute and record the stream-in speeds for each associated halo node.
         By convention, the halo nodes will be sorted by incoming/outgoing partition
         and by increasing global node number
         
        """
        self.halo_nodes_g = []   # global node number of all lattice points in this partition's halo
        self.boundary_nodes_g = [] # global node number of all lattice points on this partition's boundary 
        self.communication_list_out = [] # list of tuples containing partition, global LP and speed that must be
        # transfered from the HALO to neighboring partitions
        self.communication_list_in = [] # list of tuples containing origin partition, global LP and speed 
        # that is transfered from neighboring lattice point into boundary nodes of the partition.
        ngb_list = [] # local list of neighbors for constructing the list of NFC_Halo_Data_Organizer objects
        for l_nd in range(self.num_local_nodes):
            for spd in range(1,self.numSpd): # spd 0 is always in the local partition
                tgt_node_g = self.adjacency[l_nd,spd];
                tgt_part = self.parts[tgt_node_g];
                if tgt_part != self.rank: # the target node is in a different partition
                    ngb_list.append(tgt_part)
                    self.halo_nodes_g.append(tgt_node_g)
                    self.boundary_nodes_g.append(self.local_to_global[l_nd])
                    self.communication_list_out.append((tgt_part,tgt_node_g,spd))
                    self.communication_list_in.append((tgt_part,self.local_to_global[l_nd],self.bb_Spd[spd]))
        # when done looping through adjacency list, remove repeated elements of the list.
        self.halo_nodes_g = np.unique(self.halo_nodes_g)  # this should be sorted too.
        self.boundary_nodes_g = np.unique(self.boundary_nodes_g) 
        self.bnl_l = [] # boundary node list.  Local node number of all boundary nodes.
        for bn in self.boundary_nodes_g:
            self.bnl_l.append(self.global_to_local[bn])
        self.bnl_l = sorted(self.bnl_l[:]) # make it sorted

        """
         now that I know how many halo nodes there are, I need to assign local node numbers to them.
    
         also, need to make lists with local node numbers for interior and boundary nodes

        """
        self.num_halo_nodes = len(self.halo_nodes_g);
        self.total_nodes = self.num_local_nodes + self.num_halo_nodes;


        # assign local node numbers to the halo nodes and add to the local_to_global map
        ln = self.num_local_nodes
        for hn in self.halo_nodes_g:
            self.local_to_global[ln] = hn;
            self.global_to_local[hn] = ln;

        #print "rank %d has %d local nodes and %d halo nodes" % (self.rank, self.num_local_nodes,self.num_halo_nodes)


        """
          Let's work here to sort out the commmunication requirement.  The NFC_Part_Communicator will handle the
          actual comms, but this object needs to know the details of how to order the outgoing data and how
          to distribute the incoming data.
           
        """
        ngb_list = np.unique(ngb_list)# this will also be sorted.
        #print "rank %d has %d neighbors: %s" % (self.rank,len(ngb_list), str(ngb_list))
        self.HDO_list = [] # halo data organizer list
        for ngb in ngb_list:
            self.HDO_list.append(NFC_Halo_Data_Organizer(ngb))
       



    def load_parts(self):
        """
           read parts.lbm and get a list of lattice points that I own.
           create a global-to-local and local-to-global map of lattice points
        """

        self.parts = np.empty([self.Nx*self.Ny*self.Nz],dtype=np.int32);
        self.part_sizes = np.zeros(self.size,dtype=np.int32); # so each process knows how many
        # lattice points are local to each partition -- particularly for offseting MPI write operations.
        self.local_to_global = {}; self.global_to_local = {};
        indx = 0;  #initialize the global counter
        self.num_local_nodes = 0
        with open('parts.lbm') as parts:
            for p in parts:
                p_i = np.int32(p); # convert to np.int32
                self.part_sizes[p_i]+= 1 # increment the counter for particular partition num_local_nodes
                self.parts[indx] = p_i # store in my local array (will need to use this repeatedly)
                if p_i == self.rank: # if this lp is assigned to the current rank:
                    self.local_to_global[self.num_local_nodes] = indx; # put into local-to-global dictionary
                    self.global_to_local[indx] = self.num_local_nodes; # put in global-to-local dictionary
                    self.num_local_nodes+=1
                indx+=1 # either way increment the global counter

        

       

        
        

    def write_partition_vtk(self):
        """
          output a vtk file to visualize partition body and halo points
        """
        partition_body = np.zeros(self.Nx*self.Ny*self.Nz);
        pp = np.where(self.parts == self.rank)
        partition_body[pp] = 100 # set interior nodes for each partition
        # also signify boundary nodes and halo nodes for each rank
        partition_body[self.halo_nodes_g] = 200 # halo nodes have max value for partition
        partition_body[self.boundary_nodes_g] = 150 # boundary nodes --- intermediate to interior and halo
        vtk_filename = 'partition_map' + str(self.rank) + '.vtk'
        dims = [int(self.Nx), int(self.Ny), int(self.Nz)]
        origin = [0., 0., 0.]; 
        spacing = [1., 1., 1.];
        writeVTKpt(partition_body,'partition',vtk_filename,dims,origin,spacing)
                



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
        self.ex = [0,-1,0,0,-1,-1,-1,-1,0,0,-1,-1,-1,-1,1,0,0,1,1,1,1,0,0,1,1,1,1]
        self.ey = [0,0,-1,0,-1,1,0,0,-1,-1,-1,-1,1,1,0,1,0,1,-1,0,0,1,1,1,1,-1,-1]
        self.ez = [0,0,0,-1,0,0,-1,1,-1,1,-1,1,-1,1,0,0,1,0,0,1,-1,1,-1,1,-1,1,-1]
        self.bbSpd = [0,14,15,16,17,18,19,20,21,22,23,24,25,26,
	      1,2,3,4,5,6,7,8,9,10,11,12,13]
        self.w = [8./27.,2./27.,2./27.,2./27.,1./54.,1./54.,1./54.,1./54.,1./54.,
	       1./54.,1./216.,1./216,1./216.,1./216.,2./27.,2./27.,
	       2./27.,1./54.,1./54.,1./54.,1./54.,1./54.,
		1./54.,1./216.,1./216,1./216.,1./216.]




if __name__=="__main__":
    """
     put testing code here
    """

   
