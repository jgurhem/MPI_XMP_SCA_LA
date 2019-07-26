#!/bin/sh -e
set -e

. ~/mpi/scripts/load.sh

cd /gpfsdata/jgurhem/res/
rm -rf *.bin *.dat core.*

size=6
np=4
nrhs=1

DIR_UTILS=~/mpi/utils
DIR_EXE=~/mpi/executables
DIR_SRC=~/mpi/sources

compileSCA() {
    mpicc -Wall "$DIR_SRC"/$1 -o "$DIR_EXE"/${1%.*} "$DIR_SRC"/parse_args.c -L/gpfs1l/gpfshome/mds/staff/jgurhem/mpi/scalapack/scalapack-2.0.2/ -lscalapack -L/gpfslocal/pub/lapack/build_3.5_gnu47/lib -llapack -lblas -lifcore -Wunused-variable -O2
}

compileSCA sca_lu.c

mpirun -n $np "$DIR_EXE"/genBin $size
mpirun -n $np "$DIR_EXE"/sca_lu -s $size -A a.bin -B lu.bin
