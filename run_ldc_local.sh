#!/bin/bash

# arguments
# 
# 1 - N_divs for LDC
# 2 - lattice type [ 'D3Q15' | 'D3Q19' | 'D3Q27' ]
# 3 - partition methodology [ '1D' | '3D' | 'metis' ]
# 4 - number of partitions
# 5 - number of omp threads
# 6 - pre-process




MAT_FILE=LDC_geom_test.mat

Num_ts=50001
ts_rep_freq=50
Warmup_ts=0
plot_freq=500
Re=1000
dt=0.001
Cs=0
Restart_flag=0

# must re-process if you change:
# N_divs, partition methodology, or the number of partitions.
# if the lattice type changes, you do not have to re-do pre-processing,
# but the resulting partitions may not be the same as what would have
# been picked with new lattice type
if [ "$6" = "1" ]; then
python ./lid_driven_cavity_geom.py $1


python ./pyNFC_partition.py $MAT_FILE $2 $3 $4



else
echo "pre-processing skipped, commencing time steps"
fi

# basically, pyNFC_preprocess.py just writes params.lbm now.
python ./pyNFC_preprocess.py $MAT_FILE $2 $3 $4 \
$Num_ts $ts_rep_freq $Warmup_ts $plot_freq $Re $dt $Cs $Restart_flag

export OMP_NUM_THREADS=$5
#aprun -n $4 -d $5  ./pyNFC_run.py
mpirun -np $4 ./pyNFC_run_local.py

python ./processNFC.py 
