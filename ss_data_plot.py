#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 11:58:51 2018

@author: stu
"""

import sys
sys.path.insert(1,'.')

#import math
#import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt

# load the data
ss_filename = 'ss_data.mat' 
ss_data = scipy.io.loadmat(ss_filename)
ss_ux = ss_data['ss_ux'];
ss_uy = ss_data['ss_uy'];
ss_uz = ss_data['ss_uz'];
ss_rho = ss_data['ss_rho'];

ss_vmag = np.sqrt(ss_ux**2. + ss_uy**2. + ss_uz**2.)

im = plt.imshow(ss_uz[0,:,:],cmap='plasma')

plt.show()