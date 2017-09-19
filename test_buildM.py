# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 18:08:55 2017

@author: stu
"""

import pyLattice as pl

Nx = 1; Ny = 1; Nz = 1;

latt = pl.D3Q15Lattice(Nx,Ny,Nz)

print "Pre-orthogonalization:"
for i in range(15):
    print "row %d = "%i
    print latt.Mt[i,:]








