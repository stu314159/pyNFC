#!/bin/bash

## should take 3 arguments:
# lattice type (D3Q15 | D3Q19 | D3Q27)
# partition type  (1D | 3D | metis)
# number of partitions  (positive integer)
# number of OMP threads

# pre-process

# assume going-in that the GNU programming environment
# is selected.  May need a way to be smarter about this,
# but the goal is to change to INTEL if 'metis' is the
# selected partitioning scheme

if [ "$2" = "metis" ]; then
  module swap PrgEnv-gnu PrgEnv-intel
fi
PYTHON_EXE=/p/home/sblair/anaconda2/bin/python
$PYTHON_EXE ./pyNFC_preprocess.py $1 $2 $3

# switch back to GNU
if [ "$2" = "metis" ]; then
  module swap PrgEnv-intel PrgEnv-gnu
fi

export OMP_NUM_THREADS=$4
# invoke pyNFC
aprun -n $3 -d $4 ./pyNFC_run.py

# post-process the results
#aprun -n 10 ./pyNFC_postprocess.py
#$PYTHON_EXE ./processNFC.py
./processNFC_serial
