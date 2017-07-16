import LDH
import numpy as np
import ctypes



numSpd = 15;
myLDH15 = LDH.LBM_DataHandler(numSpd)
print "LBM_DataHandler object created with %d speeds."%myLDH15.get_numSpd()

fIn = np.array(range(0,numSpd),dtype=np.float32);

myLDH15.set_fIn(fIn)
print "fIn = ", fIn
myLDH15.multFin(2.)
print "now fIn = ", fIn


