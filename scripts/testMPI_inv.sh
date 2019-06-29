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

compileMPI mpi_inv_gj.c


coo=cooToMat
inv=cooInv
mat=binToASCII_mat_row

mpirun -n $np "$DIR_EXE"/genBin $size
mpirun -n $np "$DIR_EXE"/mpi_inv_gj $size
"$DIR_EXE"/$mat a.bin a.dat $size $size
"$DIR_EXE"/$mat inv.bin inv.dat $size $size
echo a
"$DIR_EXE"/$coo a.dat $size $size
echo inv ok
"$DIR_EXE"/$inv a.dat $size
echo inv
"$DIR_EXE"/$coo inv.dat $size $size

