import LBM_Interface as LB
import numpy as np



numSpd = 15;
myLB15 = LB.PyLBM_Interface(numSpd)
print "LBM_DataHandler object created with %d speeds."%myLB15.get_numSpd()

fIn = np.array(range(0,numSpd),dtype=np.float32);

myLB15.set_fIn(fIn)
print "fIn = ", fIn



