#genJobfile.py
"""
more-or-less automated generation of PBS jobfile
"""

import argparse

parser = argparse.ArgumentParser(prog="genJobfile.py",
                                description="PBS Jobfile generation script.")

parser.add_argument('jobfileName',type=str)
parser.add_argument('jobName',type=str)
parser.add_argument('nnodes',type=int)
parser.add_argument('ppn',type=int)
parser.add_argument('mpi_procs_per_node',type=int)
parser.add_argument('ompthreads_per_node',type=int)
parser.add_argument('runtimeNumHours',type=int)
parser.add_argument('queue',type=str)
parser.add_argument('latticeType',type=str)
parser.add_argument('partitionType',type=str)


# parse input arguments
args = parser.parse_args()

# assign to the variables

jobfileName = args.jobfileName
jobName = args.jobName
nnodes = args.nnodes
ppn = args.ppn
mpi_procs_per_node = args.mpi_procs_per_node
ompthreads_per_node=args.ompthreads_per_node
runtimeNumHours = args.runtimeNumHours
queue = args.queue
latticeType = args.latticeType
partitionType = args.partitionType

# make these additional arguments?
geom_file = 'gridChan_geom.py'
run_script = 'run_gridChan_grace.sh'
N_divs=61
pp_bool=1
dynamics=3


filesToCopy = ['FluidChannel.py', 'pyLattice.py', 'pyNFC.py', 'pyNFC_run.py',
               'pyNFC_Util.py', 'validate.py', 'vtkHelper.py', run_script,
               'pyPartition.py','pyNFC_preprocess.py','pyNFC_partition.py',
               'partition_suggestion.py','partition_compare.py',
               'LBM_Interface.so','PartitionHelper.so','processNFC.py','hdf5Helper.py',geom_file]



if runtimeNumHours < 10:
    walltime = "0%d:00:00"%runtimeNumHours
else:
    walltime = "%d:00:00"%runtimeNumHours # may be a problem if runtime > 99 hours
    
mpi_procs = mpi_procs_per_node*nnodes

jobfileName = "%s.pbs"%(jobfileName)
#--------- more-or-less fixed code below -----------------

proj_id = 'USNAM37752431'

# open the file
jf = open(jobfileName,'w')

# essential PBS directives
jf.write('#!/bin/bash \n') # the shell
#jf.write('#PBS -A %s \n'%proj_id) # project identifier
jf.write('#PBS -q %s \n'%queue) # specify queue
jf.write('#PBS -l select=%d:ncpus=%d:mpiprocs=%d:ompthreads=%d \n'% \
         (nnodes,ppn,mpi_procs_per_node,ompthreads_per_node))
jf.write('#PBS -l walltime=%s \n'%walltime)
jf.write('#PBS -l ccm=1 \n') # specify cluster compatibility mode.  Why wouldn't you?

#optional PBS directives
jf.write('#PBS -N %s \n'%jobName)
jf.write('#PBS -j oe \n')
#jf.write('#PBS -V \n')
jf.write('#PBS -S /bin/bash \n')


# Execution block
jf.write('cd /mnt/lustre/scratch/sblair/jobout \n')
jf.write("JOBID=`echo $PBS_JOBID | cut -d '.' -f 1` \n")
jf.write('if [ ! -d $JOBID ]; then \n')
jf.write('  mkdir -p $JOBID \n')

jf.write('fi \n')
jf.write('cd $JOBID \n')
# copy files 
for s in filesToCopy:
    jf.write('cp $PBS_O_WORKDIR/%s . \n'% s)

# invoke execution
jf.write('module swap PrgEnv-cray PrgEnv-gnu\n') #let's just plan on using GNU for now
#jf.write('module load costinit\n')
#jf.write('module load python\n')
#jf.write('module load numpy\n')
#jf.write('module load scipy\n')
#jf.write('module load mpi4py\n')
#jf.write('module load boost\n')
#jf.write('module load cray-hdf5\n')
jf.write('source ~/work/projects/venv/bin/activate \n');
jf.write('./%s %d %s %d %s %d %d %d\n'%(run_script,N_divs,
         latticeType,dynamics,partitionType,mpi_procs,ompthreads_per_node,pp_bool))
jf.write('deactivate\n')
jf.close()
