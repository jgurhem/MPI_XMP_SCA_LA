#@ class            = clallmds+
#@ job_name         = mpi-run
#@ total_tasks      = NCORES
#@ node             = NNODES
#@ wall_clock_limit = TIMEWALL
#@ output           = $(job_name).$(jobid).log
#@ error            = $(job_name).$(jobid).err
#@ job_type         = mpich
#@ environment      = COPY_ALL
#@ node_usage       = not_shared
#@ queue
#

. scripts/poincare/load_env.sh
rm -rf core.*
bash scripts/run.sh Poincare APP EXE NNODES NCORES DATASIZE LANG test_mpiscaxmp_DATASIZE.json
