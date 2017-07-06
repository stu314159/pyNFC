#!/bin/bash

# pre-process
python ./pyNFC_preprocess.py

# invoke pyNFC
mpirun -np 6 ./pyNFC_test.py

# post-process the results
mpirun -np 10 ./pyNFC_postprocess.py
