#!/bin/sh -e
set -e

. ~/mpi/scripts/load.sh

cd /gpfsdata/jgurhem/res/
rm -rf *.bin *.dat core.*

size=8
np=4
nrhs=1

DIR_UTILS=~/mpi/utils
DIR_EXE=~/mpi/executables
DIR_SRC=~/mpi/sources

compile() {
    xmpcc -DCOO_OUT -Wall "$DIR_SRC"/$1 -o "$DIR_EXE"/${1%.*}
}

compile xmp_sls_lu.c


mpirun -n $np "$DIR_EXE"/genBin $size
mpirun -n $np "$DIR_EXE"/xmp_sls_lu $size

check_results -op slsg -one-file -s 1 -b $size -A a.dat -V v.dat -R r.dat -ff coo -print

