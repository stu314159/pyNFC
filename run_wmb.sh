#!/bin/bash

# arguments
# 
# 1 - N_divs for wall mounted brick 
# 2 - lattice type [ 'D3Q15' | 'D3Q19' | 'D3Q27' ]
# 3 - partition methodology [ '1D' | '3D' | 'metis' ]
# 4 - number of partitions
# 5 - number of omp threads
# 6 - pre-process





if [ "$6" = "1" ]; then
python ./wmb_geom.py $1

# saves mat file named ChanCavityTest.mat
MAT_FILE=wall_mounted_brick.mat

Num_ts=201
ts_rep_freq=50
Warmup_ts=0
plot_freq=20
Re=25
dt=0.01
Cs=0
Restart_flag=0

python ./pyNFC_preprocess.py $MAT_FILE $2 $3 $4 \
$Num_ts $ts_rep_freq $Warmup_ts $plot_freq $Re $dt $Cs $Restart_flag

else
echo "pre-processing skipped, commencing time steps"
fi

export OMP_NUM_THREADS=$5
mpirun -np $4 ./pyNFC_run.py

python ./processNFC.py 


