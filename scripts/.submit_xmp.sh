#@ class            = clallmds+
#@ job_name         = out-run
#@ total_tasks      = PROCS
#@ node             = NODES
#@ wall_clock_limit = 05:00:00
#@ output           = $(job_name).$(jobid).log
#@ error            = $(job_name).$(jobid).err
#@ job_type         = mpich
#@ environment      = COPY_ALL
#@ queue
#

. ~/mpi/scripts/load.sh

path=$(pwd)

cd /gpfsdata/jgurhem/res/
rm -f *.bin

DIR_EXE=/gpfshome/mds/staff/jgurhem/mpi/executables

mpirun -n $LOADL_TOTAL_TASKS "$DIR_EXE"/genBin SIZE
bash "$path/run.sh" blockGauss SIZE "$DIR_EXE"/xmp_sls_g NODES PROCS XMP
bash "$path/run.sh" blockGaussJordan SIZE "$DIR_EXE"/xmp_sls_gj NODES PROCS XMP
bash "$path/run.sh" blockLU SIZE "$DIR_EXE"/xmp_lu NODES PROCS XMP
bash "$path/run.sh" blockLUsolveLS SIZE "$DIR_EXE"/xmp_sls_lu NODES PROCS XMP

