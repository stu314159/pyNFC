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
# 8 - restart
# 9 - time average


# saves mat file named ChanCavityTest.mat
MAT_FILE=wall_mounted_brick.mat

Num_ts=1001
ts_rep_freq=500
Warmup_ts=0
plot_freq=100
Re=25
dt=0.01
Cs=0
Restart_flag=$8
TimeAvg_flag=$9

if [ "$7" = "1" ]; then
aprun -n 1 ./wmb_geom.py $1

aprun -n 1 ./pyNFC_partition.py $MAT_FILE $2 $4 $5

else
echo "pre-processing skipped, commencing time steps"
fi

aprun -n 1 ./pyNFC_preprocess.py $MAT_FILE $2 $3 $4 $5 \
$Num_ts $ts_rep_freq $Warmup_ts $plot_freq $Re $dt $Cs $Restart_flag \
$TimeAvg_flag

export OMP_NUM_THREADS=$6

aprun -n $5 -d $6  ./pyNFC_run.py

aprun -n 1 ./processNFC.py 


