#!/bin/sh -e
set -e

DIR_UTILS=~/mpi/utils
DIR_EXE=~/mpi/executables
DIR_SRC=~/mpi/sources

. ~/mpi/scripts/load.sh

compileGCC() {
    gcc -o "$DIR_EXE"/${1%.*} "$DIR_UTILS"/$1 -std=c99
}

compileSCA() {
    mpicc -Wall "$DIR_SRC"/$1 -o "$DIR_EXE"/${1%.*} -L/gpfs1l/gpfshome/mds/staff/jgurhem/mpi/scalapack/scalapack-2.0.2/ -lscalapack -L/gpfslocal/pub/lapack/build_3.5_gnu47/lib -llapack -lblas -lifcore -Wunused-variable -O2
}

compileMPI() {
    mpicc -Wall "$DIR_SRC"/$1 -o "$DIR_EXE"/${1%.*}
}

compileGCC cooToMat.c
compileGCC binToASCII_mat_col.c
compileGCC binToASCII_mat_row.c

compileSCA lu.c
compileSCA sls.c
compileSCA inv.c
compileSCA genBin.c

compileMPI gauss.c
compileMPI gaussJordan.c
