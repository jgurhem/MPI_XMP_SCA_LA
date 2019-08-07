#!/bin/sh -e
set -e

. poincare/load_env.sh

cd /gpfsdata/jgurhem/res/
rm -rf *.bin *.dat core.*

size=8
np=4
nrhs=1

DIR_UTILS=~/mpi/utils
DIR_EXE=~/mpi/executables
DIR_SRC=~/mpi/sources

compile() {
    xmpcc -DCOO_OUT -Wall "$DIR_SRC"/$1 "$DIR_SRC"/mpiio_dmat.c "$DIR_SRC"/parse_args.c -o "$DIR_EXE"/${1%.*}
}

compile xmp_lu.c

mpirun -n $np "$DIR_EXE"/genBin $size
mpirun -n $np "$DIR_EXE"/xmp_lu -s $size -A a.bin -B lu.bin

check_results -op blu -one-file -s 1 -b $size -A a.dat -B lu.dat -ff coo -print
check_results -op blu -one-file -s 1 -b $size -A a.bin -B lu.bin -ff binR -print
