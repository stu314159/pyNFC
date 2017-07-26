#!/bin/bash

## should take 4 arguments:
# 1 lattice type (D3Q15 | D3Q19 | D3Q27)
# 2 partition type  (1D | 3D | metis)
# 3 number of partitions  (positive integer)
# 4 number of OMP threads

# pre-process
python ./pyNFC_preprocess.py $1 $2 $3

export OMP_NUM_THREADS=$4
# invoke pyNFC
mpirun -n $3 ./pyNFC_run.py

# post-process the results
python ./processNFC.py
