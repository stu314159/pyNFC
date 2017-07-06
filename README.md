# pyNFC
Python implementation of Numeric Fluid Channel

This is under development and, frankly, is quite crude.  The goal of this project is to 
serve as a prototype for planned improvements to the Numeric Fluid Channel code residing
in nfc2.  

Current status: appears to be working properly (albiet slowly) for all lattice
types with the metis partitioning.  Scalability has not been formally evaluated
but appears to be quite strong (near linear scaling for single problem size)

Working now to improve performance on a per-process basis while maintaining 
correctness.
