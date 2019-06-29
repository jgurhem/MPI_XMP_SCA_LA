#!/bin/sh -e
set -e

DIR_UTILS=~/mpi/utils
DIR_EXE=~/mpi/executables
DIR_SRC=~/mpi/sources

rm -fv "$DIR_EXE"/*

. ~/mpi/scripts/load.sh

compileGCC() {
    gcc -o "$DIR_EXE"/${1%.*} "$DIR_UTILS"/$1 -std=c99
}

compileSCA() {
    mpicc -Wall "$DIR_SRC"/$1 -o "$DIR_EXE"/${1%.*} -L/gpfs1l/gpfshome/mds/staff/jgurhem/mpi/scalapack/scalapack-2.0.2/ -lscalapack -L/gpfslocal/pub/lapack/build_3.5_gnu47/lib -llapack -lblas -lifcore -Wunused-variable -O2
}

compileMPI() {
    mpicc -Wall "$DIR_SRC"/$1 "$DIR_SRC"/mpiio_dmat.c -o "$DIR_EXE"/${1%.*}
}

compileXMP() {
    xmpcc -Wall "$DIR_SRC"/$1 "$DIR_SRC"/mpiio_dmat.c -o "$DIR_EXE"/${1%.*}
    rm -f ${1%.*}.o
}

compileGCC cooToMat.c
compileGCC binToASCII_mat_col.c
compileGCC binToASCII_mat_row.c
compileGCC cooInv.c
compileGCC cooLU.c

compileSCA sca_lu.c
compileSCA sca_sls_lu.c
compileSCA sca_inv_lu.c
compileSCA genBin.c

compileMPI mpi_sls_g.c
compileMPI mpi_sls_gj.c
compileMPI mpi_sls_lu.c
compileMPI mpi_lu.c
compileMPI mpi_inv_gj.c

compileXMP xmp_sls_g.c
compileXMP xmp_sls_gj.c
compileXMP xmp_sls_lu.c
compileXMP xmp_lu.c
