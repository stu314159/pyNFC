#!/bin/bash

## should take 3 arguments:
# lattice type (D3Q15 | D3Q19 | D3Q27)
# partition type  (1D | 3D | metis)
# number of partitions  (positive integer)

# pre-process
python ./pyNFC_preprocess.py $1 $2 $3

# invoke pyNFC
mpirun -np $3 ./pyNFC_test.py

# post-process the results
python ./processNFC.py
