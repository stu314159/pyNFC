# pyNFC
Python implementation of Numeric Fluid Channel

The goal of this project is to serve as a prototype for planned improvements to the 
Numeric Fluid Channel code residing in nfc2.  

Current status: working and tested thoroughly for D3Q15 lattice type with metis
partitioning.  Uses Boost.Python to achieve significant performance improvments
over Python, Numpy, and mpi4py alone.  Scalability is under evaluation and appears
to be quite good.  Single node performance (with 32 MPI processes) is slightly better
than NFCpp on a single node with OpenMP acceleration.

Working now to more formally validate D3Q19 and D3Q27 lattice structures as well as
evaluate pro/con of metis partitions.
