#@ class            = clallmds+
#@ job_name         = out-run
#@ total_tasks      = PROCS
#@ node             = NODES
#@ wall_clock_limit = 20:00:00
#@ output           = $(job_name).$(jobid).log
#@ error            = $(job_name).$(jobid).err
#@ job_type         = mpich
#@ environment      = COPY_ALL 
#@ queue
#

. ~/mpi/scripts/load.sh
#source /gpfs1l/opt/Intel/itac/8.1.0.024/bin/itacvars.sh
cd /gpfsdata/jgurhem/res/
rm -f *.bin

DIR_EXE=/gpfshome/mds/staff/jgurhem/mpi/executables

execRun(){

t=( $(mpirun -n $LOADL_TOTAL_TASKS $3 $2) )
echo ${t[*]}
success="true"

if [ ${#t[*]} -ne 3 ]
then
	t="error"
	success="false"
fi

echo "Poincare;$LOADL_TOTAL_TASKS;NODES;$1;SCA;$2;$LOADL_TOTAL_TASKS;$(date);-1;0;0;${t[0]};0;$success;;Compilation -O2" >> ~/results.csv
echo "Poincare;$LOADL_TOTAL_TASKS;NODES;$1;SCAL;$2;$LOADL_TOTAL_TASKS;$(date);-1;0;0;${t[1]};0;$success;;Compilation -O2" >> ~/results.csv
echo "Poincare;$LOADL_TOTAL_TASKS;NODES;$1;SCALS;$2;$LOADL_TOTAL_TASKS;$(date);-1;0;0;${t[2]};0;$success;;Compilation -O2" >> ~/results.csv

}

mpirun -n $LOADL_TOTAL_TASKS "$DIR_EXE"/genBin SIZE
execRun ScalapackLU SIZE "$DIR_EXE"/sca_lu
execRun ScalapackLU_sls SIZE "$DIR_EXE"/sca_sls_lu
execRun ScalapackLU_inv SIZE "$DIR_EXE"/sca_inv_lu

