#!/bin/bash

# pre-process
python ./pyNFC_preprocess.py

# invoke pyNFC
aprun -n 36 ./pyNFC_test.py

# post-process the results
aprun -n 10 ./pyNFC_postprocess.py
