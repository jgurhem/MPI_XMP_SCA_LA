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

compileMPI() {
    mpicc -Wall "$DIR_SRC"/$1 -o "$DIR_EXE"/${1%.*}
}

compileMPI gauss.c


coo=cooToMat
mat=binToASCII_mat_col

mpirun -n $np "$DIR_EXE"/genBin $size
mpirun -n $np "$DIR_EXE"/gauss $size
"$DIR_EXE"/$mat a.bin a.dat $size $size
"$DIR_EXE"/$mat b.bin b.dat $size $nrhs
"$DIR_EXE"/$mat r.bin r.dat $size $nrhs
echo a
"$DIR_EXE"/$coo a.dat $size $size
echo b
"$DIR_EXE"/$coo b.dat $size $nrhs
echo r
"$DIR_EXE"/$coo r.dat $size $nrhs

