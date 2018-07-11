#!/bin/bash

# arguments
#
# 1 - N_divs for wall mounted brick 
# 2 - lattice type [ 'D3Q15' | 'D3Q19' | 'D3Q27' ]
# 3 - dynamics [ 1 = LBGK | 2 = RBGK | 3 = MRT]
# 4 - partition methodology [ '1D' | '3D' | 'metis' ]
# 5 - number of partitions
# 6 - number of omp threads
# 7 - pre-process <-- for now, assume we will always pre-process
# 8 - restart
# 9 - time average


# saves mat file named ChanCavityTest.mat
MAT_FILE=gridChan.mat

Num_ts=10001
ts_rep_freq=1000
Warmup_ts=0
plot_freq=5000
Re=15000
dt=0.000005
Cs=20
Restart_flag=$8
TimeAvg_flag=$9

# pre-processing job will be done with only
# one process
pre_job=$(qsub gridChan_Onyx_pre.pbs)

# dependent "run" job will invoke more resources
run_job=$(qsub -W depend=afterok:$pre_job gridChan_Onyx_run.pbs)

# dependent "post" job will, again, only require a single process
post_job=$(qsub -w depende=afterok:$run_job gridChan_Onyx_post.pbs)
