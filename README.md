# pyNFC
Python implementation of Numeric Fluid Channel.  

Distributed parallel with mpi4py.  Calculations mainly done by C++-based LBM library coupled with Python code through Boost.Python.  OpenMP provides some acceleration and may be used in conjunction with MPI.  

Dependencies: 
(Required)
1. NumPy
2. SciPy
3. MPI4Py
4. Boost.Python
5. h5py 

(Optional)
6. pymetis

Allows for LBM simulation of fluid flows.  Bulk dynamics include: LBGK, RBGK, and MRT.  Lattice Options: D3Q15, D3Q19, and D3Q27 (as of now, no MRT for D3Q27).  Single parameter turbulence model implemented for LBGK and RBGK.  Regularized velocity and pressure boundary conditions implemented for all lattices on the west (Z-min) and east (Z-max) boundaries respectively.  Boundary conditions also included for a simple moving boundary (only in the z-direction) and, of course, solid boundaries.

Code not updated for a while.  The main issue is getting the build/operating environment set on HPC systems.  Currently not able to get the mpi, boost, and python environments right for the LBM_Interface build to work.  tLBM is my principal LBM project currently and only depends on MPI and HDF5.
