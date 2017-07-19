#!/bin/bash

## should take 3 arguments:
# lattice type (D3Q15 | D3Q19 | D3Q27)
# partition type  (1D | 3D | metis)
# number of partitions  (positive integer)

# pre-process
PYTHON_EXE=/p/home/sblair/anaconda2/bin/python
$PYTHON_EXE ./pyNFC_preprocess.py $1 $2 $3

# invoke pyNFC
aprun -n $3 ./pyNFC_test.py

# post-process the results
#aprun -n 10 ./pyNFC_postprocess.py
$PYTHON_EXE ./processNFC.py
