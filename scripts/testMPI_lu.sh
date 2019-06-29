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

compileMPI() {
    mpicc -Wall "$DIR_SRC"/$1 "$DIR_SRC"/mpiio_dmat.c -o "$DIR_EXE"/${1%.*}
}

compileMPI mpi_lu.c

mat=binToASCII_mat_col

mpirun -n $np "$DIR_EXE"/genBin $size
mpirun -n $np "$DIR_EXE"/mpi_lu $size
"$DIR_EXE"/$mat a.bin a.dat $size $size
"$DIR_EXE"/$mat lu.bin lu.dat $size $size

check_results -op blu -one-file -s 1 -b $size -A a.dat -B lu.dat -ff coo -print
