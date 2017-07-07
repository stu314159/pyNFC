#pyNFC.py
"""
Implementation file for pyNFC library module -- under development

"""

import numpy as np
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTKpt
from mpi4py import MPI

import pyLattice as pl
from pyNFC_Util import NFC_Halo_Data_Organizer

class NFC_LBM_partition(object):
    """
    each partition has:
         
    """
    def __init__(self,rank,size,comm,Nx,Ny,Nz,rho_lbm,u_bc,omega,Cs,lattice_type='D3Q15',):
        """
          rank - MPI rank for this partition
          size - MPI size for MPI COMM WORLD
          comm - MPI communicator
          Nx, Ny, Nz so the partition has info about the overall
                     lattice structure
          lattice_type - ['D3Q15' | 'D3Q19' | 'D3Q27']

          rho_lbm - scaled density for outlet boundary condition
          u_bc - scaled velocity for inlet boundary condition
          omega - relaxation constant for LBM collisions
          Cs - parameter for turbulence model
        """
        self.rank = rank; self.size = size; self.comm = comm # MPI variabes
        self.Nx = Nx; self.Ny = Ny; self.Nz = Nz; # all partitions need to know the global domain structure

        # LBM simulation parameters
        self.rho_lbm = rho_lbm; self.u_bc = u_bc; self.omega = omega; self.Cs = Cs
        
        
        if lattice_type == 'D3Q15':
            self.lattice = pl.D3Q15Lattice(self.Nx, self.Ny, self.Nz)
        elif lattice_type == 'D3Q19':
            self.lattice = pl.D3Q19Lattice(self.Nx, self.Ny, self.Nz)
        else:
            self.lattice = pl.D3Q27Lattice(self.Nx, self.Ny, self.Nz)

        self.numSpd = self.lattice.get_numSpd()
        #print "process %d of %d constructed %s lattice " % (rank,size,lattice_type)
        self.ex = np.array(self.lattice.get_ex(),dtype=np.int32);
        self.ey = np.array(self.lattice.get_ey(),dtype=np.int32);
        self.ez = np.array(self.lattice.get_ez(),dtype=np.int32);
        self.bb_Spd = np.array(self.lattice.get_bbSpd(),dtype=np.int32);
        self.w = np.array(self.lattice.get_w(),dtype=np.float32);

        # these functions must be called in sequence to ensure the appropriate
        # member data has been constructed
        self.load_parts()
        self.gen_adjacency() # adjacency list using global node numbers
        self.get_halo_nodes() # halo nodes are all global node numbers - local node number lists also produced
        self.get_interior_nodes() # local node numbers of all interior nodes produced
        self.allocate_data_arrays() # all global data arrays that will be required for the LBM simulation
        self.initialize_node_lists() # snl, inl and onl lists from pre-processed data

        # halo nodes are now incorporated into local node lists and data arrays.
        # convert the adjacency matrix so that it is local
        self.convert_adjacency_to_local() 

        self.initialize_data_arrays() # fEven, fOdd - in future, include restart data load.
        self.out_requests = [MPI.REQUEST_NULL for i in range(self.num_ngb)]
        self.in_requests = [MPI.REQUEST_NULL for i in range(self.num_ngb)]
        self.statuses = [MPI.Status() for i in range(self.num_ngb)]
        

        # mpi file writing variables
        self.vtk_dump_num = 0;
        self.vtk_ux_stub = 'ux'; self.vtk_uy_stub = 'uy'; self.vtk_uz_stub = 'uz'
        self.vtk_rho_stub = 'density'

        self.vtk_suffix = '.b_dat'


       
        
        #self.report_statistics() # say something about yourself
        # comment out when you are done visualizing the partitions 
        #self.write_partition_vtk() # visualize each partition interior, boundary and halo


    def take_LBM_timestep(self,isEven):
        """
          carry out the LBM process for a time step.  Orchestrate processing of all
          lattice points and communicating data between MPI partitions.
        """
        
        # process boundary lattice points
        #tst_rank = 0;
#        if self.rank == tst_rank:
#            print "rank %d processing boundary nodes" % (tst_rank)   
        self.process_lattice_points(isEven,self.bnl_l)

        # extract halo data
#        if self.rank == tst_rank:
#            print "rank %d extracting halo data" % (tst_rank)
        self.extract_halo_data(isEven)

        # initiate communication of halo data
#        if self.rank == tst_rank:
#            print "rank %d initiate send/recv of halo data" % (tst_rank)
        for ngb in range(self.num_ngb):
            ngb_rnk = self.ngb_list[ngb]
            self.out_requests[ngb] = self.comm.Isend([self.HDO_out_dict[ngb_rnk].buffer,
                                            self.HDO_out_dict[ngb_rnk].buff_len,
                                            MPI.FLOAT],ngb_rnk,self.rank)
            self.in_requests[ngb] = self.comm.Irecv([self.HDO_in_dict[ngb_rnk].buffer,
                                           self.HDO_in_dict[ngb_rnk].buff_len,
                                           MPI.FLOAT],ngb_rnk,MPI.ANY_TAG)
        
        

        # process interior lattice points
        #print "rank %d processing %d nodes on the interior"%(self.rank, len(self.int_l))
#        if self.rank == tst_rank:
#            print "rank %d processing interior nodes" % (tst_rank)
        self.process_lattice_points(isEven,self.int_l)

        # be sure MPI communication is done
#        if self.rank == tst_rank:
#            print "rank %d waiting for MPI coms" % (tst_rank)
        MPI.Request.Waitall(self.in_requests,self.statuses)


        # load incoming data to appropriate array
#        if self.rank == tst_rank:
#            print "rank %d inserting boundary data" % (tst_rank)
        self.insert_boundary_data(isEven)

        # done.
        

    def process_lattice_points(self,isEven,lp_list):
        """
          carry out the LBM process for a list of lattice points
          isEven - boolean to indicate if this is an even time step or odd time step
          lp_list - list of lattice points to be processed (e.g. interior lattice points or boundary lattice points)

        """

        # point fIn and fOut to the right data arrays
        if isEven:
            fIn = self.fEven; fOut = self.fOdd
        else:
            fIn = self.fOdd; fOut = self.fEven

        # initially, implement exactly as you would for C/C++
        # goal: be sure to get the correct answer; worry about performance later

        
        for lp in lp_list:
        

            # get microscopic densities
            f = fIn[lp,:]
            ndType = 0
            # get node type
            if self.inl[lp] == 1:
                ndType = 2
            elif self.onl[lp] == 1:
                ndType = 3
            elif self.snl[lp] == 1:
                ndType = 1
            
                
             
            # process lattice point and get outlet value
            
            f_o = self.lattice.compute_fOut(f,ndType,self.omega,self.Cs,self.u_bc,self.rho_lbm)

            # stream to outlet value
            self.stream(fOut,f_o,lp)

           
    def stream(self,fOut,f,lp):
        """
            stream collided particle density distributions to neighbor lattice points
        """

        for spd in range(self.numSpd):
            tgtNd = self.adjacency[lp,spd]
            fOut[tgtNd,spd] = f[spd]


    def initialize_node_lists(self):
        """
         load pre-processor data into inl, onl and snl lists
        """
        
        inl_filename = "inl.lbm"; 
        
        inl_f = open(inl_filename,'r');
        numINL = int(inl_f.readline());
        for i in range(numINL):
            gIN = int(inl_f.readline());
            if gIN in self.global_to_local:
               lIN = self.global_to_local[gIN]
               self.inl[lIN] = 1
        inl_f.close()

        onl_filename = "onl.lbm";
        onl_f = open(onl_filename,'r');
        numONL = int(onl_f.readline());
        for i in range(numONL):
            gOUT = int(onl_f.readline());
            if gOUT in self.global_to_local:
                lOUT = self.global_to_local[gOUT]
                self.onl[lOUT] = 1
        onl_f.close()

        snl_filename = "snl.lbm";
        snl_f = open(snl_filename,'r');
        numSNL = int(snl_f.readline());
        for i in range(numSNL):
            gSNL = int(snl_f.readline());
            if gSNL in self.global_to_local:
                lSNL = self.global_to_local[gSNL]
                self.snl[lSNL] = 1
        snl_f.close()

    def write_node_sorting(self):
        """
           write a file "ordering.b_dat" that will contain the global 
           node numbers of each local node (in order)

           this will facilitate post-processing

        """
        # create numpy array with global node numbers of local nodes
        node_roster = np.empty([self.num_local_nodes],dtype=np.int32)

        # poplulate the node roster from self.num_local_nodes + self.local_to_global
        for nd in range(self.num_local_nodes):
            node_roster[nd] = self.local_to_global[nd]


        # file mode
        amode = MPI.MODE_WRONLY | MPI.MODE_CREATE
        
        # file name
        file_name = 'ordering.b_dat'

        #
        fh = MPI.File.Open(self.comm,file_name,amode)
        offset = self.offset_int*np.dtype(np.int32).itemsize
        fh.Write_at_all(offset,node_roster) 
        fh.Close()


    def write_data(self,isEven):
        """
          write your partition of data to the MPI data file

        """

#        if self.rank == 0:
#            print "computing local macroscopic data"

        ux, uy, uz, rho = self.compute_local_data(isEven);

#        if self.rank == 0:
#            print "writing data"

        # self.offset_bytes is the number of bytes offset

        # file mode
        amode = MPI.MODE_WRONLY | MPI.MODE_CREATE

        # create file names
        ux_fn = self.vtk_ux_stub + str(self.vtk_dump_num) + self.vtk_suffix
        uy_fn = self.vtk_uy_stub + str(self.vtk_dump_num) + self.vtk_suffix
        uz_fn = self.vtk_uz_stub + str(self.vtk_dump_num) + self.vtk_suffix
        rho_fn = self.vtk_rho_stub + str(self.vtk_dump_num) + self.vtk_suffix
        
        # write ux
        fh = MPI.File.Open(self.comm,ux_fn,amode)
        fh.Write_at_all(self.offset_bytes,ux)
        fh.Close()

        # write uy
        fh = MPI.File.Open(self.comm,uy_fn,amode)
        fh.Write_at_all(self.offset_bytes,uy)
        fh.Close()

        # write uz
        fh = MPI.File.Open(self.comm,uz_fn,amode)
        fh.Write_at_all(self.offset_bytes,uz)
        fh.Close()

        # write rho
        fh = MPI.File.Open(self.comm,rho_fn,amode)
        fh.Write_at_all(self.offset_bytes,rho)
        fh.Close()
        

        # when done writing, increment vtk dump number
        self.vtk_dump_num += 1


    def compute_local_data(self,isEven):
        """
          compute ux, uy, uz and rho for output at requested data dump intervals
        """
        if isEven:
            f = self.fEven;
        else:
            f = self.fOdd;

        # prepare arrays to hold the data
        ux = np.zeros([self.num_local_nodes],dtype=np.float32)
        uy = np.zeros([self.num_local_nodes],dtype=np.float32)
        uz = np.zeros([self.num_local_nodes],dtype=np.float32)
        rho = np.zeros([self.num_local_nodes],dtype=np.float32)
    

        for lp in self.lnl_l: # this sucks.  basically requires that 
        # self.lnl_l goes from 0 to self.num_local_nodes...
            for spd in range(self.numSpd):
                rho[lp]+=f[lp,spd];
                ux[lp]+=self.ex[spd]*f[lp,spd];
                uy[lp]+=self.ey[spd]*f[lp,spd];
                uz[lp]+=self.ez[spd]*f[lp,spd];
            ux[lp]/=rho[lp]; uy[lp]/=rho[lp]; uz[lp]/=rho[lp]
        ux[np.where(self.snl[:self.num_local_nodes]==1)] = 0.;
        uy[np.where(self.snl[:self.num_local_nodes]==1)]= 0.;
        uz[np.where(self.snl[:self.num_local_nodes]==1)] = 0.;
        uz[np.where(self.inl[:self.num_local_nodes]==1)] = self.u_bc # set boundary condition
        
                
            

        return ux, uy, uz, rho
        


    def report_statistics(self):
        """
          gather and report a collection of interesting (?) statistics:
              - number of local nodes
              - number of interior nodes and boundary nodes
              - number of ngb partitions
              - number of f_alphas to communicate in total

        """

        num_local_nodes = self.num_local_nodes
        num_interior_nodes = len(self.int_l)
        num_boundary_nodes = len(self.bnl_l)
        num_nbg = len(self.HDO_out_dict)
        num_dat = 0
        for key in self.HDO_out_dict:
            num_dat+=self.HDO_out_dict[key].count_data_members()

        print "rank %d has %d local nodes, %d interior nodes, %d boundary nodes, %d neighbors and %d data elements to commmunicate"%(self.rank, num_local_nodes, num_interior_nodes, num_boundary_nodes, num_nbg, num_dat)
    

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
                
                

    def convert_adjacency_to_local(self):
        """
          convert the elements of the adjacency matrix into local node numbers

        """
        for nd in range(self.num_local_nodes):
            for spd in range(self.numSpd):
                lInd = self.global_to_local[self.adjacency[nd,spd]]
                self.adjacency[nd,spd] = lInd


    def get_interior_nodes(self):
        """
         generate a sorted list of all nodes not on the boundary or halo.
         assumes that self.get_halo_nodes(self) has already been called by the
         constructor so that halo and boundary node lists have already been formed

        """
        all_nodes = np.arange(self.total_nodes,dtype=np.int32)
        not_interior = np.union1d(np.array(self.bnl_l,dtype=np.int32),
                                  np.array(self.hnl_l,dtype=np.int32))

        self.int_l = np.setxor1d(all_nodes,not_interior) # self.int_l = interior node list in local node numbers


    def get_halo_nodes(self):
        """
         Identify global node numbers that must be incorporated onto the halo.
         Identify the stream-out speeds associated with each halo node and the partition
         each halo node belongs in.
         Also compute and record the stream-in speeds for each associateself.halo_nodes_gd halo node.
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

        """
        assign local node numbers to the halo nodes and add to the local_to_global map
        assume halo nodes are given local node numbers in order of increasing global node
        number starting after all local node numbers.  This may not be optimal, so be 
        prepared to change it.
        """
        ln = self.num_local_nodes
        for hn in self.halo_nodes_g:
            self.local_to_global[ln] = hn;
            self.global_to_local[hn] = ln;
            ln+=1 # increment the local node counter so halo nodes are assigned unique numbers

        #print "rank %d has %d local nodes and %d halo nodes" % (self.rank, self.num_local_nodes,self.num_halo_nodes)
        self.lnl_l = range(self.num_local_nodes); # "local node list - local (node numbers) "again, be prepared to change this.
        # make a list of the local node numbers of all halo nodes
        self.hnl_l = []  # halo node list.  (in case I use a different scheme for locally numbering halo nodes)
        for hn in self.halo_nodes_g:
            self.hnl_l.append(self.global_to_local[hn])

        """
          Let's work here to sort out the commmunication requirement.  The NFC_Part_Communicator will handle the
          actual comms, but this object needs to know the details of how to order the outgoing data and how
          to distribute the incoming data.
           
        """
        ngb_list = np.unique(ngb_list)# this will also be sorted.
        #print "rank %d has %d neighbors: %s" % (self.rank,len(ngb_list), str(ngb_list))
        self.HDO_out_dict = {} # halo data organizer dictionary for outgoing
        self.HDO_in_dict = {} # halo data organizer dictionary for incoming data
        for ngb in ngb_list:
            self.HDO_out_dict[ngb] = NFC_Halo_Data_Organizer(ngb)
            self.HDO_in_dict[ngb] = NFC_Halo_Data_Organizer(ngb)
        #iterate through communication list to put the communication items into the correct HDO
        for c in self.communication_list_out:
            self.HDO_out_dict[c[0]].insert(c[1],c[2])
        
        for c in self.communication_list_in:
            self.HDO_in_dict[c[0]].insert(c[1],c[2])

        # now construct global and local node lists required for data book-keeping
        for ngb in ngb_list:
            self.HDO_out_dict[ngb].make_lists()
            self.HDO_out_dict[ngb].make_lists_local(self.global_to_local)
            self.HDO_out_dict[ngb].allocate_buffer() # buffer is _.buffer
            self.HDO_in_dict[ngb].make_lists()
            self.HDO_in_dict[ngb].make_lists_local(self.global_to_local)
            self.HDO_in_dict[ngb].allocate_buffer()


        self.ngb_list = ngb_list[:] # make this a data member
        self.num_ngb = len(self.ngb_list)
            
       
    def extract_halo_data(self,isEven):
        """
           function to extract data from halo nodes and 
           place into data buffers for HDO objects

        """
        if isEven:
            fOut = self.fOdd
        else:
            fOut = self.fEven

        for ngb in self.ngb_list:
            self.HDO_out_dict[ngb].extract_halo_data(fOut)


    def insert_boundary_data(self,isEven):
        """
          function to insert data passed from neighboring partitions
          into the boundary points where they belong
        """

        if isEven:
            f = self.fOdd
        else:
            f = self.fEven

        for ngb in self.ngb_list:
            self.HDO_in_dict[ngb].insert_boundary_data(f)


    def load_parts(self):
        """
           read parts.lbm and get a list of lattice points that I own.
           create a global-to-local and local-to-global map of lattice points

           these maps are later updated to include halo nodes associated with each partition

           lastly - get a count of lattice points associated with each partition so that 
           the local rank can calculate the necessary offset for file outputs
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
                    self.global_to_local[indx] = self.num_local_nodes; # put in global-to-localbml dictionary
                    self.num_local_nodes+=1
                indx+=1 # either way increment the global counter


        # save the cumsum of all partitions with rank lower than self.rank
        # to use in offsetting MPI write operations.
        self.offset_int = np.sum(self.part_sizes[0:self.rank]) # this should exclude the current rank.
        offset_check = np.sum(self.part_sizes[0:self.rank+1])
        assert (offset_check - self.offset_int) == self.num_local_nodes
        
        self.offset_bytes = self.offset_int * np.dtype(np.float32).itemsize;

        
        

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
                


    def allocate_data_arrays(self):
        """
         allocate arrays for LBM simulation
        """
        # some thought/testing should be done regarding the shape of this data array.
        self.fEven = np.empty([self.total_nodes , self.numSpd],dtype=np.float32)
        self.fOdd = np.empty_like(self.fEven)
        self.snl = np.zeros([self.total_nodes],dtype=np.int32);
        self.inl = np.zeros([self.total_nodes],dtype=np.int32);
        self.onl = np.zeros([self.total_nodes],dtype=np.int32);

    def initialize_data_arrays(self):
        """
         set initial density distributions.  Currently implemented to
         initialize everything to the outlet pressure boundary condition with
         zero velocity

        """
        
        for idx in range(self.total_nodes):
            self.fEven[idx,:] = self.rho_lbm * self.w
            self.fOdd[idx,:] = self.rho_lbm * self.w





if __name__=="__main__":
    """
     put testing code here
    """

   
