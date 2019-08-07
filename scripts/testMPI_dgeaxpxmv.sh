#!/bin/sh
#set -e

. poincare/load_env.sh

cd /gpfsdata/jgurhem/res/
rm -rf *.bin *.dat core.*

size=10
np=4
nrhs=1

DIR_UTILS=~/mpi/utils
DIR_EXE=~/mpi/executables
DIR_SRC=~/mpi/sources

compileMPI() {
    mpicc -Wall "$DIR_SRC"/$1 "$DIR_SRC"/mpiio_dmat.c "$DIR_SRC"/parse_args.c -o "$DIR_EXE"/${1%.*}
}

compileMPI mpi_dgeaxpxmv.c

mpirun -n $np "$DIR_EXE"/genBin $size
mpirun -n $np "$DIR_EXE"/mpi_dgeaxpxmv -s $size -A a.bin -V b.bin -R r.bin

check_results -op dgeaxpxmv -one-file -s 1 -b $size -A a.bin -V b.bin -R r.bin -ff binR -print

