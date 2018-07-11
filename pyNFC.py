#pyNFC.py
"""
Implementation file for pyNFC library module -- under development

"""

import numpy as np
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTKpt
from mpi4py import MPI

import pyLattice as pl
from pyNFC_Util import NFC_Halo_Data_Organizer
import LBM_Interface as LB
import h5py
import scipy.io

class NFC_LBM_partition(object):
    """
    each partition has:
         
    """
    def __init__(self,rank,size,comm,Nx,Ny,Nz,rho_lbm,u_bc,dynamics,omega,Cs,lattice_type='D3Q15',):
        """
          rank - MPI rank for this partition
          size - MPI size for MPI COMM WORLD
          comm - MPI communicator
          Nx, Ny, Nz so the partition has info about the overall
                     lattice structure
          lattice_type - ['D3Q15' | 'D3Q19' | 'D3Q27']

          rho_lbm - scaled density for outlet boundary condition
          u_bc - scaled velocity for inlet boundary condition
          dynamics - 1 = LBGK | 2 = RBGK | 3 = MRT
          omega - relaxation constant for LBM collisions
          Cs - parameter for turbulence model
        """
        self.rank = rank; self.size = size; self.comm = comm # MPI variabes
        self.Nx = Nx; self.Ny = Ny; self.Nz = Nz; # all partitions need to know the global domain structure

        # LBM simulation parameters

        self.rho_lbm = rho_lbm; self.u_bc = u_bc; self.omega = omega; self.Cs = Cs
        self.dynamics = dynamics;
        self.timeAvg = False

        
        
        if lattice_type == 'D3Q15':
            self.lattice = pl.D3Q15Lattice(self.Nx, self.Ny, self.Nz)
            
        elif lattice_type == 'D3Q19':
            self.lattice = pl.D3Q19Lattice(self.Nx, self.Ny, self.Nz)
        else:
            self.lattice = pl.D3Q27Lattice(self.Nx, self.Ny, self.Nz)

              
                
        self.numSpd = self.lattice.get_numSpd()
        self.myLB = LB.PyLBM_Interface(self.numSpd) # boost interface
        self.myLB.set_Ubc(self.u_bc)
        self.myLB.set_rhoBC(self.rho_lbm)
        self.myLB.set_omega(self.omega)
        self.myLB.set_dynamics(self.dynamics)
        self.myLB.set_Cs(self.Cs)
        self.myLB.set_MPIcomm(self.comm)
        self.myLB.set_timeAvg(self.timeAvg)
        
        # if dynamics == 3, construct lattice.omegaMRT and pass its pointer to myLB
        if self.dynamics == 3:
            self.lattice.set_omega(omega)
            self.lattice.constructOmegaMRT(self.omega);
            self.myLB.set_omegaMRT(self.lattice.omegaMRT); # pass the MRT operator pointer to myLB
        
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
        self.initialize_node_lists() # ndType, ssNds, and ssNd_list 

        # halo nodes are now incorporated into local node lists and data arrays.
        # convert the adjacency matrix so that it is local
        self.convert_adjacency_to_local() 

        self.initialize_data_arrays() # fEven, fOdd - in future, include restart data load.
        self.out_requests = [MPI.REQUEST_NULL for i in range(self.num_ngb)]
        self.in_requests = [MPI.REQUEST_NULL for i in range(self.num_ngb)]
        self.statuses = [MPI.Status() for i in range(self.num_ngb)]
        
        # pass pointers of node lists to myLB object
#    
        self.myLB.set_ndT(self.ndT) #replace use of inl, onl and snl
        #self.myLB.set_ssNds(self.ssNds) # assign ss Node list
        self.myLB.set_ssNds(self.ss_nd_array) # assign ss Node lis
        self.myLB.set_adjacency(self.adjacency)
        self.myLB.set_fEven(self.fEven)
        self.myLB.set_fOdd(self.fOdd)
        self.myLB.set_boundaryNL(self.bnl_l)
        self.myLB.set_bnlSZ(int(self.num_bn))
        self.myLB.set_inlSZ(int(self.num_in))
        self.myLB.set_interiorNL(self.int_l)
        self.myLB.set_totalNodes(int(self.total_nodes))
        
        

        # mpi file writing variables
        self.vtk_dump_num = 0;
        self.vtk_ux_stub = 'ux'; self.vtk_uy_stub = 'uy'; self.vtk_uz_stub = 'uz'
        self.vtk_rho_stub = 'density'

        self.vtk_suffix = '.b_dat'


    def initialize_timeAvg(self):
        """
        if time-averaging is requested, set boolean specifying time averaging
        and allocate data arrays to hold time-average data.
        
        Pointers to the data array need to be passed to the PyLBM_Interface object
        """
        self.timeAvg = True
                     
        self.uAvg = np.zeros([self.num_local_nodes],dtype=np.float32);
        self.vAvg = np.zeros([self.num_local_nodes],dtype=np.float32);
        self.wAvg = np.zeros([self.num_local_nodes],dtype=np.float32);
        self.rhoAvg = np.zeros([self.num_local_nodes],dtype=np.float32);
        
        # pass pointers to PyLBM_Interface object
        self.myLB.set_uAvg(self.uAvg);
        self.myLB.set_vAvg(self.vAvg);
        self.myLB.set_wAvg(self.wAvg);
        self.myLB.set_rhoAvg(self.rhoAvg);
        self.myLB.set_timeAvg(True);
        
    def write_timeAvg(self):
        """
        when the simulation is done:
        a) divide all time average data by the number of time steps; and
        b) write the local data to disk in appropriately-named data files.
        
        """
        
                
        # self.uAvg/=float(NumTs);
        # self.vAvg/=float(NumTs);
        # self.wAvg/=float(NumTs);
        # self.rhoAvg/=float(NumTs);
        
        # self.offset_bytes is the number of bytes offset
        
        # file mode
        amode = MPI.MODE_WRONLY | MPI.MODE_CREATE

        # create file names
        ux_fn = 'uAvg.b_dat'
        uy_fn = 'vAvg.b_dat'
        uz_fn = 'wAvg.b_dat'
        rho_fn = 'rhoAvg.b_dat'
        
                       
        # write uAvg
        fh = MPI.File.Open(self.comm,ux_fn,amode)
        fh.Write_at_all(self.offset_bytes,self.uAvg[:self.num_local_nodes])
        fh.Close()

        # write vAvg
        fh = MPI.File.Open(self.comm,uy_fn,amode)
        fh.Write_at_all(self.offset_bytes,self.vAvg[:self.num_local_nodes])
        fh.Close()

        # write wAvg
        fh = MPI.File.Open(self.comm,uz_fn,amode)
        fh.Write_at_all(self.offset_bytes,self.wAvg[:self.num_local_nodes])
        fh.Close()

        # write rhoAvg
        fh = MPI.File.Open(self.comm,rho_fn,amode)
        fh.Write_at_all(self.offset_bytes,self.rhoAvg[:self.num_local_nodes])
        fh.Close()
        
    def take_LBM_timestep(self,isEven):
        """
          carry out the LBM process for a time step.  Orchestrate processing of all
          lattice points and communicating data between MPI partitions.
        """
        
        # process boundary lattice points
 
        self.myLB.process_nodeList(isEven,0);

        # extract halo data

        self.myLB.extract_halo_data(isEven)

        # initiate communication of halo data

        for ngb in range(self.num_ngb):
            ngb_rnk = self.ngb_list[ngb]
            self.out_requests[ngb] = self.comm.Isend([self.HDO_out_dict[ngb_rnk].buffer,
                                            self.HDO_out_dict[ngb_rnk].buff_len,
                                            MPI.FLOAT],ngb_rnk,self.rank)
            self.in_requests[ngb] = self.comm.Irecv([self.HDO_in_dict[ngb_rnk].buffer,
                                           self.HDO_in_dict[ngb_rnk].buff_len,
                                           MPI.FLOAT],ngb_rnk,MPI.ANY_TAG)
        
        

        # process interior lattice points
  
        self.myLB.process_nodeList(isEven,1)

        # be sure MPI communication is done

        MPI.Request.Waitall(self.in_requests,self.statuses)


        # load incoming data to appropriate array

        self.myLB.insert_boundary_data(isEven)

        # done.
        

   

           
    def stream(self,fOut,f,lp):
        """
            stream collided particle density distributions to neighbor lattice points
        """

        #for spd in range(self.numSpd):
        #    tgtNd = self.adjacency[lp,spd]
        #    fOut[tgtNd,spd] = f[spd]
        tgtNds = self.adjacency[lp,:]
        fOut[tgtNds,range(self.numSpd)]=f[:]


    def initialize_node_lists(self):
        """
         load pre-processor data into ndType lists
        """
        

        
        ndType_filename = "ndType.lbm";
        ndl_f = open(ndType_filename,'r');
        numNt = int(self.Nx*self.Ny*self.Nz)
        for i in range(numNt): # iterate through all members of the ndType file
            gNt = int(ndl_f.readline()); # get the node type for the global node number
            if i in self.global_to_local: # check if this global node number is in my partition
                lNd = self.global_to_local[i] #if it is, get the local node number
                self.ndT[lNd] = gNt # set the ndT list to the indicated number
        ndl_f.close()
        
        # need to read "ssNds.lbm" to determine which of my nodes are ssNds.
        
        self.part_ss_sizes = np.zeros(self.size,dtype=np.int32); # so
        # each process can know how many subspace nodes other partitions have
        
        ssNds_filename = "ssNds.lbm"
        self.ss_nd_list = []
        ss_f = open(ssNds_filename)
        for line in ss_f.readlines():
            gNt = int(line); # get the global node number
            self.part_ss_sizes[self.parts[gNt]]+=1
            if gNt in self.global_to_local: # if this global node number is in this partition
                lNd = self.global_to_local[gNt] # get the local partition index
                if (lNd <= self.num_local_nodes):
                    self.ssNds[lNd] = 1 # set the subspace node list value to 1
                    self.ss_nd_list.append(lNd)
                    
        ss_f.close()
        #print "rank %d has %d subspace nodes"%(self.rank,len(self.ss_nd_list))
        self.ss_nd_list.sort() # sort the list for convenience.
        self.ss_nd_array = np.array(self.ss_nd_list,dtype=np.int32) # make this into a numpy array for Boost.
        self.ss_offset_int = np.sum(self.part_ss_sizes[0:self.rank]);
        #print "rank %d has offset int equal to %d"%(self.rank, self.ss_offset_int);
        
        # if rank == 0, write the part_ss_sizes to disk.<<<------****
        if self.rank==0:
            filename = 'part_ss_sizes.mat'
            part_dict = {}
            part_dict['part_ss_sizes']=list(self.part_ss_sizes[:])
            scipy.io.savemat(filename,part_dict)
            
            
    def write_ss_node_sorting(self):
        """
          write a file "ss_ordering.b_dat" that will contain the global
          node numbers of each subspace node (in order) <-- mpi file
          
          
        """
        ss_node_roster = np.empty([len(self.ss_nd_list)],dtype=np.int32);
        for nd in range(len(self.ss_nd_list)):
            ss_node_roster[nd]=self.local_to_global[self.ss_nd_list[nd]]
            
        amode = MPI.MODE_WRONLY | MPI.MODE_CREATE
        file_name = 'ss_ordering.b_dat'
        fh = MPI.File.Open(self.comm,file_name,amode)
        offset = self.ss_offset_int*np.dtype(np.int32).itemsize;
        fh.Write_at_all(offset,ss_node_roster)
        fh.Close();
            
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

    def write_subspace_data(self):
        """
        at the end of the simulation, write all accumulated subspace data to 
        binary data files for subsequent processing
        """
        amode = MPI.MODE_WRONLY | MPI.MODE_CREATE
        
        # create file names
        ss_ux_fn = "ss_ux.b_dat"
        ss_uy_fn = "ss_uy.b_dat"
        ss_uz_fn = "ss_uz.b_dat"
        ss_rho_fn = "ss_rho.b_dat"
        
        self.ss_offset_bytes = self.num_ts*self.ss_offset_int*(np.dtype(np.float32).itemsize)
        
        # write ss_ux
        fh = MPI.File.Open(self.comm,ss_ux_fn,amode)
        fh.Write_at_all(self.ss_offset_bytes,self.ssNd_ux);
        fh.Close()
        
        # write ss_uy
        fh = MPI.File.Open(self.comm,ss_uy_fn,amode)
        fh.Write_at_all(self.ss_offset_bytes,self.ssNd_uy);
        fh.Close()
        
        # write ss_uz
        fh = MPI.File.Open(self.comm,ss_uz_fn,amode)
        fh.Write_at_all(self.ss_offset_bytes,self.ssNd_uz);
        fh.Close()
        
        # write ss_rho
        fh = MPI.File.Open(self.comm,ss_rho_fn,amode)
        fh.Write_at_all(self.ss_offset_bytes,self.ssNd_rho);
        fh.Close()

    def write_data(self,isEven):
        """
          write your partition of data to the MPI data file

        """

        ux, uy, uz, rho = self.compute_local_data(isEven);

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
       
        # prepare arrays to hold the data
        ux = np.zeros([self.num_local_nodes],dtype=np.float32)
        uy = np.zeros([self.num_local_nodes],dtype=np.float32)
        uz = np.zeros([self.num_local_nodes],dtype=np.float32)
        rho = np.zeros([self.num_local_nodes],dtype=np.float32)
    
        self.myLB.set_ux(ux);
        self.myLB.set_uy(uy);
        self.myLB.set_uz(uz);
        self.myLB.set_rho(rho);

        self.myLB.compute_local_data(isEven);

       
        ux[np.where(self.ndT[:self.num_local_nodes]==1)] = 0.;
        uy[np.where(self.ndT[:self.num_local_nodes]==1)]= 0.;
        uz[np.where(self.ndT[:self.num_local_nodes]==1)] = 0.;
        uz[np.where(self.ndT[:self.num_local_nodes]==2)] = self.u_bc # set boundary condition
        
                
            

        return ux, uy, uz, rho
        
    def record_subspace_data(self,ts):
        """
        compute ux, uy, uz, and rho for subspace data set and store in the appropriate
        copy of the subspace data arrays.
        """
        self.myLB.compute_subspace_data(ts)

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
        self.num_in = self.int_l.size; # number of interior nodes

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
        self.num_bn = len(self.bnl_l)
        self.bnl_l = np.array(self.bnl_l,dtype=np.int32)

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
        
        # register and initialize neighbors in PyLBM_Interface object
        for ngb in self.ngb_list:
            numData = self.HDO_out_dict[ngb].count_data_members()
            self.myLB.registerNeighbor(int(ngb),int(numData))
            # get pointers to incoming halo data
            self.myLB.getHaloInPointers(self.HDO_in_dict[ngb].lnn_array,
                                        self.HDO_in_dict[ngb].spd_array,
                                        self.HDO_in_dict[ngb].buffer,int(ngb))
            # get pointers to outgoing halo data
            self.myLB.getHaloOutPointers(self.HDO_out_dict[ngb].lnn_array,
                                        self.HDO_out_dict[ngb].spd_array,
                                        self.HDO_out_dict[ngb].buffer,int(ngb))
            
       
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
        self.ndT = np.zeros([self.total_nodes],dtype=np.int32);
        self.ssNds = np.zeros([self.total_nodes],dtype=np.int32);
        

    def allocate_subspace_data_arrays(self,num_ts):
        """
        Allocate data arrays for storing subset data for this partition. over 
        all time steps.  At the end of the simulation, all of the subspace
        data will be written to a single binary data file using MPI's 
        API.  This data file then will need to be post processed so that it can
        be visualized and or further processed for turbulence data.
        
        """
        #print "rank %d allocating subspace data arrays for %d timesteps"%(self.rank,num_ts)
        self.ssNd_ux = np.zeros([num_ts,len(self.ss_nd_list)],dtype=np.float32)
        self.ssNd_uy = np.zeros_like(self.ssNd_ux)
        self.ssNd_uz = np.zeros_like(self.ssNd_ux)
        self.ssNd_rho = np.zeros_like(self.ssNd_ux)
        
        # assign pointers through Boost Interface
        self.myLB.set_ss_ux(self.ssNd_ux)
        self.myLB.set_ss_uy(self.ssNd_uy)
        self.myLB.set_ss_uz(self.ssNd_uz)
        self.myLB.set_ss_rho(self.ssNd_rho)
        
        self.myLB.set_num_ssNds(len(self.ss_nd_list));
        
        self.num_ts = num_ts;
        
        
        
    def initialize_data_arrays(self):
        """
         set initial density distributions.  Currently implemented to
         initialize everything to the outlet pressure boundary condition with
         zero velocity

        """
        
        for idx in range(self.total_nodes):
            self.fEven[idx,:] = self.rho_lbm * self.w
            self.fOdd[idx,:] = self.rho_lbm * self.w

    def load_restart_data(self):
        """
         load velocity and density data from restart.h5 and use data
         to initialize fEven and fOdd arrays.  Note that this is called
         from pyNFC_run.py and it occurs after construction and initialization
         of the fEven and fOdd data arrays.
        """
        
        # # allocate numpy arrays to store the local node data
        # ux_l = np.zeros((self.total_nodes,1),dtype=np.float);
        # uy_l = np.zeros_like(ux_l);
        # uz_l = np.zeros_like(ux_l);
        # rho_l = np.zeros_like(ux_l);
        
        # open the restart.h5 file for reading
        f = h5py.File('restart.h5','r')
        ux_d = f['velocity/x'];
        uy_d = f['velocity/y'];
        uz_d = f['velocity/z'];
        rho_d = f['density/rho'];
        
        for nd in range(self.total_nodes):
            gNd = self.local_to_global[nd]
            # ux_l[nd] = ux_d[gNd];
            # uy_l[nd] = uy_d[gNd];
            # uz_l[nd] = uz_d[gNd];
            # rho_l[nd] = rho_d[gNd];
            u = [ux_d[gNd],uy_d[gNd],uz_d[gNd]]
            rho = rho_d[gNd];
            self.fEven[nd,:] = self.lattice.compute_equilibrium([],rho,u);
            self.fOdd[nd,:] = self.lattice.compute_equilibrium([],rho,u);
                      
        # when done, close h5py
        f.close()
        
        # convert ux, uy, uz, rho data into density distribution values
        # for fEven and fOdd --> load into the arrays.



if __name__=="__main__":
    """
     put testing code here
    """

   
