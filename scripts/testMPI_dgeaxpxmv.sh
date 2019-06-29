#!/bin/sh
#set -e

. ~/mpi/scripts/load.sh

cd /gpfsdata/jgurhem/res/
rm -rf *.bin *.dat core.*

size=10
np=4
nrhs=1

DIR_UTILS=~/mpi/utils
DIR_EXE=~/mpi/executables
DIR_SRC=~/mpi/sources

compileMPI() {
    mpicc -Wall "$DIR_SRC"/$1 "$DIR_SRC"/mpiio_dmat.c -o "$DIR_EXE"/${1%.*}
}

compileMPI mpi_dgeaxpxmv.c


coo=cooToMat
mat=binToASCII_mat_row

echo generating data
mpirun -n $np "$DIR_EXE"/genBin $size
echo dgaxpxmv
mpirun -n $np "$DIR_EXE"/mpi_dgeaxpxmv $size
echo converting a
"$DIR_EXE"/$mat a.bin a.dat $size $size
echo converting b
"$DIR_EXE"/$mat b.bin b.dat $size $nrhs
echo converting r
"$DIR_EXE"/$mat r.bin r.dat $size $nrhs
echo a
"$DIR_EXE"/$coo a.dat $size $size
echo b
"$DIR_EXE"/$coo b.dat $size $nrhs
echo r
"$DIR_EXE"/$coo r.dat $size $nrhs

