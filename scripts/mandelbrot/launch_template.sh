rm -rf core.*

. scripts/mandelbrot/load_env.sh

CORES_PER_NODES=$(( NCORES / NNODES ))

{
for node in $(seq 1 NNODES)
do
  for i in $(seq 1 $CORES_PER_NODES)
  do
    echo "mandelbrot-${node}"
  done
done
} > machinefile.txt

bash scripts/run.sh Mandelbrot APP "-machinefile machinefile.txt EXE" NNODES NCORES DATASIZE LANG test_mpiscaxmp_DATASIZE.json
