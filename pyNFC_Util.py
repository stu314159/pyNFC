# pyNFC_Util.py
"""
definition module for some pyNFC utility classes


"""

import numpy as np

class NFC_Halo_Data_Organizer(object):
    """
     collect and organize how halo data will be organized prior to communication
    """
    def __init__(self,ngb_rank):
        """
          constructor
        """
        self.ngb_rank = ngb_rank
        self.halo_data = {}

    def get_ngb_rank(self):
        return self.ngb_rank

    def insert(self,gnn,spd):
        """
         data insertion function:
         gnn - global node number of outgoing/incoming halo data
         spd - speed of outgoing/incoming halo data
        """
        self.halo_data.setdefault(gnn,[]).append(spd) # more pythonic
        

    def make_lists(self): # only need be done once
        """
          once all data is inserted, this function will create
          two lists: the gnn for outgoing data to this rank; and
                     the vector of speeds for the outgoing data.
        """
        self.gnn_list = []
        self.spd_list = []
        sorted_keys = sorted(self.halo_data.keys()); # sort the gnn keys
        for k in sorted_keys:
            values = self.halo_data[k]; values = sorted(values) # sort the spds for each gnn
            for v in values:
                self.gnn_list.append(k)
                self.spd_list.append(v)

              

    def make_lists_local(self,global_to_local):
        """
         make and store local node number version of lists (but keep gnn ordering) -- be careful
        """
        self.lnn_list = []
        for g in self.gnn_list: # preserve order?
            self.lnn_list.append(global_to_local[g]) #preserve order?

    def count_data_members(self):
        """
          report the number of items to be sent from this partition
        """
        return len(self.gnn_list)

    def allocate_buffer(self):
        """
          construct buffer for data in/out
        """
        self.buff_len = len(self.gnn_list)
        self.buffer = np.empty([self.buff_len],dtype=np.float32)


    def extract_halo_data(self,fOut):
        """
          extract required data from fOut and place into buffer
        """
        for d in range(len(self.lnn_list)):
            ln = self.lnn_list[d]; spd = self.spd_list[d];
            self.buffer[d] = fOut[ln,spd]

    def insert_boundary_data(self,f):
        """
           insert stream-in data into the appropriate boundary node/speed
        """
        for d in range(len(self.lnn_list)):
            ln = self.lnn_list[d]; spd = self.spd_list[d];
            f[ln,spd] = self.buffer[d]



#class NFC_Part_Communicator(object):
#    """
#     class designed to handle the communication tasks for an NFC_LBM_partition
#    """
#    def __init__(self,rank,size,comm,comm_list_gnn, comm_list_spd):
#        """
#          rank - which MPI process
#          size - number of MPI processes
#          comm - MPI communicator
#          comm_list_gnn - list of global node numbers
#          
#
#          Each partition needs to know:
#            a) the list of MPI ranks it needs to exchange halo data with;
#            b) for each exchange pair the # of data elements to be exchanged;
#            
#            Of course each process also has to manage the local/global node numbers
#            for which each data element is bound as well as the associated speed.
#
#            (for both incoming and outgoing data) as well as what order it will arrive in.
#            by convention: 
#
#                a) each partition will send data to neighboring partitions in order of 
#                   increasing global lattice point number; and
#                b) for lattice points receiving multiple speeds from the SAME neighbor partition,
#                   the data will be provided by increasing spd number.
#
#            how this is accomplished, and which classes accomplish this, is the design question
#            yet to be answered.
#          
#        """
#        self.rank = rank; self.size = size; self.comm = comm;lattice_type='D3Q15',
#        self.comm_list = comm_list;
#        self.lattice = lattice;
#        



