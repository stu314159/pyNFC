#!/bin/bash

# arguments
# 
# 1 - N_divs for wall mounted brick 
# 2 - lattice type [ 'D3Q15' | 'D3Q19' | 'D3Q27' ]
# 3 - dynamics [ 1 = LBGK | 2 = RBGK | 3 = MRT]
# 4 - partition methodology [ '1D' | '3D' | 'metis' ]
# 5 - number of partitions
# 6 - number of omp threads
# 7 - pre-process [1 = yes, 0 = no]
# 8 - Restart [1= yes, 0 = no]
# 9 - Time Avg [1 = yes, 0 = no]
# 10 - Subspace Data [1 = yes, 0 = no]



MAT_FILE=ssTest.mat


Num_ts=21
ts_rep_freq=2
Warmup_ts=0
plot_freq=5

Re=25
dt=0.01
Cs=0
Restart_flag=$8
TimeAvg_flag=$9
SubspaceData_flag=${10}


# must re-process if you change:
# N_divs, partition methodology, or the number of partitions.
# if the lattice type changes, you do not have to re-do pre-processing,
# but the resulting partitions may not be the same as what would have
# been picked with new lattice type
if [ "$7" = "1" ]; then
python ./ssTest_geom.py $1


python ./pyNFC_partition.py $MAT_FILE $2 $4 $5



else
echo "pre-processing skipped, commencing time steps"
fi

# basically, pyNFC_preprocess.py just writes params.lbm now.
python ./pyNFC_preprocess.py $MAT_FILE $2 $3 $4 $5 \
$Num_ts $ts_rep_freq $Warmup_ts $plot_freq $Re $dt $Cs $Restart_flag \
$TimeAvg_flag $SubspaceData_flag

export OMPI_MCA_mpi_warn_on_fork=0
export OMP_NUM_THREADS=$6
##aprun -n $5 -d $6  ./pyNFC_run.py
mpirun -np $5 ./pyNFC_run_numba.py

python ./processNFC.py 

#if [ "${10}" = "1" ]; then
#echo "processing subspace data"
#./process_subspace_data.py
#fi
