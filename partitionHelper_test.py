# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 08:31:14 2017

@author: stu
"""

import PartitionHelper as PH


numSpd = 15;
Nx = 5;
Ny = 5;
Nz = 5;
myPH = PH.PartitionHelper(Nx,Ny,Nz,numSpd);

adjDict = {};
myPH.setAdjacency(adjDict)
print "adjDict after setting adjacency:"
print adjDict