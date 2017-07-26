#!/bin/bash

# arguments
# 
# 1 - N_divs for wall mounted brick 
# 2 - lattice type [ 'D3Q15' | 'D3Q19' | 'D3Q27' ]
# 3 - partition methodology [ '1D' | '3D' | 'metis' ]
# 4 - number of partitions
# 5 - number of omp threads
# 6 - pre-process


# saves mat file named ChanCavityTest.mat
MAT_FILE=ChanCavityTest.mat

Num_ts=5001
ts_rep_freq=100
Warmup_ts=0
plot_freq=500
Re=100
dt=0.005
Cs=0
Restart_flag=0


# must re-process if you change:
# N_divs, partition methodology, or the number of partitions.
# if the lattice type changes, you do not have to re-do pre-processing,
# but the resulting partitions may not be the same as what would have
# been picked with new lattice type
if [ "$6" = "1" ]; then
python ./channel_cavity_geom.py $1

if [ "$3" = "metis" ]; then
  module swap PrgEnv-gnu PrgEnv-intel
fi
python ./pyNFC_partition.py $MAT_FILE $2 $3 $4


if [ "$3" = "metis" ]; then
  module swap PrgEnv-intel PrgEnv-gnu
fi

else
echo "pre-processing skipped, commencing time steps"
fi

# basically, pyNFC_preprocess.py just writes params.lbm now.
python ./pyNFC_preprocess.py $MAT_FILE $2 $3 $4 \
$Num_ts $ts_rep_freq $Warmup_ts $plot_freq $Re $dt $Cs $Restart_flag

export OMP_NUM_THREADS=$5
aprun -n $4 -d $5  ./pyNFC_run.py

#python ./processNFC.py 
./processNFC_serial

