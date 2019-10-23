#!/bin/sh -e
set -e

DIR_ROOT=$PWD

DIR_UTILS=${DIR_ROOT}/utils
DIR_EXE=${DIR_ROOT}/executables
DIR_SRC=${DIR_ROOT}/sources

mkdir -p "$DIR_EXE"
rm -fv "$DIR_EXE"/*

compileGCC() {
    gcc -o "$DIR_EXE"/${1%.*} "$DIR_UTILS"/$1 -std=c99
}

compileSCA() {
    mpicc -Wall "$DIR_SRC"/$1 -o "$DIR_EXE"/${1%.*} "$DIR_SRC"/parse_args.c -L/gpfs1l/gpfshome/mds/staff/jgurhem/mpi/scalapack/scalapack-2.0.2/ -lscalapack -L/gpfslocal/pub/lapack/build_3.5_gnu47/lib -llapack -lblas -lifcore -Wunused-variable -O2
}

compileMPI() {
    mpicc -Wall "$DIR_SRC"/$1 "$DIR_SRC"/mpiio_dmat.c "$DIR_SRC"/parse_args.c -o "$DIR_EXE"/${1%.*}
}

compileXMP() {
    xmpcc -Wall "$DIR_SRC"/$1 "$DIR_SRC"/mpiio_dmat.c "$DIR_SRC"/parse_args.c -o "$DIR_EXE"/${1%.*}
    rm -f ${1%.*}.o
}

compileGCC binToASCII_mat_col.c

compileMPI genBin.c
compileMPI mpi_sls_g.c
compileMPI mpi_sls_gj.c
compileMPI mpi_sls_lu.c
compileMPI mpi_lu.c
compileMPI mpi_inv_gj.c

