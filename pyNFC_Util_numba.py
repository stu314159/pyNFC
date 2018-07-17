# pyNFC_Util_numba.py
"""
definition module for some pyNFC utility classes with coordination with numba


"""

import numpy as np
#import numba
from numba import cuda

@cuda.jit('void(float32[:,:],int32[:],int32[:],float32[:],int32)')
def extract_halo(fOut,spd_array,lnn_array,buff,N):
    """
    kernel to get halo data from fOut and load into buff
    """
    tx = cuda.threadIdx.x;
    bx = cuda.blockIdx.x;
    bw = cuda.blockDim.x;
    tid = tx + bx*bw;
    
    if (tid < N):
        lnn = lnn_array[tid];
        spd = spd_array[tid];
        buff[tid] = fOut[lnn,spd];

@cuda.jit('void(float32[:,:],int32[:],int32[:],float32[:],int32)')
def insert_halo(fIn,spd_array,lnn_array,buff,N):
    """
    kernel to insert halo data from buff to fIn
    """
    tx = cuda.threadIdx.x;
    bx = cuda.blockIdx.x;
    bw = cuda.blockDim.x;
    tid = tx + bx*bw;
    
    if (tid < N):
        lnn = lnn_array[tid];
        spd = spd_array[tid];
        fIn[lnn,spd] = buff[tid];
    


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
        self.TPB = 128;

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
        self.d_spd_array = cuda.to_device(self.spd_array);

              

    def make_lists_local(self,global_to_local):
        """
         make and store local node number version of lists (but keep gnn ordering) -- be careful
        """
        self.lnn_list = []
        for g in self.gnn_list: # preserve order?
            self.lnn_list.append(global_to_local[g]) #preserve order?
            
        self.lnn_array = np.array(self.lnn_list,dtype=np.int32)
        self.d_lnn_array = cuda.to_device(self.lnn_array);

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
        self.d_buffer = cuda.to_device(self.buffer);
        self.nData = np.int32(len(self.lnn_list));


    def extract_halo_data(self,fOut):
        """
          extract required data from fOut and place into buffer
          
          for this version, we expect fOut to be device pointer
        """
        
        
        
        num_blocks = np.ceil(float(self.nData)/float(self.TPB))
        griddim = int(num_blocks);
        blockdim = self.TPB
        
        extract_halo[griddim,blockdim](fOut,self.d_spd_array, self.d_lnn_array, 
                    self.d_buffer, self.nData);
        
        self.buffer = self.d_buffer.copy_to_host();
        
        #for d in range(len(self.lnn_list)):
        #    ln = self.lnn_list[d]; spd = self.spd_list[d];
        #    self.buffer[d] = fOut[ln,spd]
        
        # give this a shot
        #self.buffer[:] = fOut[self.lnn_list[:],self.spd_list[:]]

    def insert_boundary_data(self,f):
        """
           insert stream-in data into the appropriate boundary node/speed
           -
           here we assume that f is a device array
        """
        # for incoming data, copy the buffere to the device. (must be done every time)
        self.d_buffer = cuda.to_device(self.buffer);
        # d_spd_array and d_lnn_array stay on device (and are unchanged)
        num_blocks = np.ceil(float(self.nData)/float(self.TPB))
        griddim = int(num_blocks);
        blockdim = self.TPB
        
        insert_halo[griddim,blockdim](f,self.d_spd_array,self.d_lnn_array,
                   self.d_buffer,self.nData);
        
        
        for d in range(len(self.lnn_list)):
            ln = self.lnn_list[d]; spd = self.spd_list[d];
            f[ln,spd] = self.buffer[d]
        
        # give this a shot
        #f[self.lnn_list[:],self.spd_list[:]]







