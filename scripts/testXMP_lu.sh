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
    xmpcc -DCOO_OUT -Wall "$DIR_SRC"/$1 "$DIR_SRC"/mpiio_dmat.c -o "$DIR_EXE"/${1%.*}
}

compile xmp_lu.c


coo=cooToMat
lu=cooLU
mat=binToASCII_mat_col

mpirun -n $np "$DIR_EXE"/genBin $size
mpirun -n $np "$DIR_EXE"/xmp_lu $size

check_results -op blu -one-file -s 1 -b $size -A a.dat -B lu.dat -ff coo -print
