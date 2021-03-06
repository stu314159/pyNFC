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
                
        self.spd_array = np.array(self.spd_list,dtype=np.int32)

              

    def make_lists_local(self,global_to_local):
        """
         make and store local node number version of lists (but keep gnn ordering) -- be careful
        """
        self.lnn_list = []
        for g in self.gnn_list: # preserve order?
            self.lnn_list.append(global_to_local[g]) #preserve order?
            
        self.lnn_array = np.array(self.lnn_list,dtype=np.int32)

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
        
        # give this a shot
        #self.buffer[:] = fOut[self.lnn_list[:],self.spd_list[:]]

    def insert_boundary_data(self,f):
        """
           insert stream-in data into the appropriate boundary node/speed
        """
        for d in range(len(self.lnn_list)):
            ln = self.lnn_list[d]; spd = self.spd_list[d];
            f[ln,spd] = self.buffer[d]
        
        # give this a shot
        #f[self.lnn_list[:],self.spd_list[:]]







