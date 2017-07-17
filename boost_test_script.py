import LBM_Interface as LB
import numpy as np



numSpd = 15;
ndType = 3;
myLB15 = LB.PyLBM_Interface(numSpd)
print "LBM_DataHandler object created with %d speeds."%myLB15.get_numSpd()

myLB15.set_ndType(ndType);

print "node type set to %d."%myLB15.get_ndType()

fIn = np.array(range(0,numSpd),dtype=np.float32);

myLB15.set_fIn(fIn)
print "fIn = ", fIn



