#!/bin/bash

# arguments
# 
# 1 - N_divs for wall mounted brick 
# 2 - lattice type [ 'D3Q15' | 'D3Q19' | 'D3Q27' ]
# 3 - dynamics [ 1 = LBGK | 2 = RBGK | 3 = MRT]
# 4 - partition methodology [ '1D' | '3D' | 'metis' ]
# 5 - number of partitions
# 6 - number of omp threads
# 7 - pre-process




MAT_FILE=wall_mounted_brick.mat

Num_ts=5001
ts_rep_freq=500
Warmup_ts=0
plot_freq=500
Re=25
dt=0.01
Cs=0
Restart_flag=0

# must re-process if you change:
# N_divs, partition methodology, or the number of partitions.
# if the lattice type changes, you do not have to re-do pre-processing,
# but the resulting partitions may not be the same as what would have
# been picked with new lattice type
if [ "$7" = "1" ]; then
python ./wmb_geom.py $1


python ./pyNFC_partition.py $MAT_FILE $2 $4 $5



else
echo "pre-processing skipped, commencing time steps"
fi

# basically, pyNFC_preprocess.py just writes params.lbm now.
python ./pyNFC_preprocess.py $MAT_FILE $2 $3 $4 $5 \
$Num_ts $ts_rep_freq $Warmup_ts $plot_freq $Re $dt $Cs $Restart_flag

export OMP_NUM_THREADS=$6
#aprun -n $5 -d $6  ./pyNFC_run.py
mpirun -np $5 ./pyNFC_run_local.py

python ./processNFC.py 
