#!/bin/sh -e
set -e

. poincare/load_env.sh

cd /gpfsdata/jgurhem/res/
rm -rf *.bin *.dat core.*

size=15
np=4
nrhs=1

DIR_UTILS=~/mpi/utils
DIR_EXE=~/mpi/executables
DIR_SRC=~/mpi/sources

compileMPI() {
    mpicc -Wall "$DIR_SRC"/$1 "$DIR_SRC"/mpiio_dmat.c -o "$DIR_EXE"/${1%.*}
}

compileMPI mpi_sls_lu.c

mat=binToASCII_mat_col

mpirun -n $np "$DIR_EXE"/genBin $size
mpirun -n $np "$DIR_EXE"/mpi_sls_lu $size | sort
#mpirun -n $np "$DIR_EXE"/mpi_sls_lu $size
"$DIR_EXE"/$mat a.bin a.dat $size $size
"$DIR_EXE"/$mat lu.bin lu.dat $size $size
"$DIR_EXE"/$mat r.bin r.dat $size $nrhs
"$DIR_EXE"/$mat b.bin b.dat $size $nrhs

check_results -op slsg -one-file -s 1 -b $size -A a.dat -V b.dat -R r.dat -ff coo -print
