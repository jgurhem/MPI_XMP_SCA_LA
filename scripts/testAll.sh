#!/bin/sh -e

. ~/mpi/scripts/load.sh

cd /gpfsdata/jgurhem/res/

size=9
np=6
nrhs=1

DIR_EXE=~/mpi/executables

coo=cooToMat
mat=binToASCII_mat_col

mpirun -n $np "$DIR_EXE"/genBin $size

mpirun -n $np "$DIR_EXE"/lu $size
"$DIR_EXE"/$mat a.bin a.dat $size $size
"$DIR_EXE"/$mat lu.bin lu.dat $size $size
echo a
"$DIR_EXE"/$coo a.dat $size $size
echo lu
"$DIR_EXE"/$coo lu.dat $size $size

echo ""
echo ""
mpirun -n $np "$DIR_EXE"/sls $size
"$DIR_EXE"/$mat a.bin a.dat $size $size
"$DIR_EXE"/$mat lu.bin lu.dat $size $size
"$DIR_EXE"/$mat b.bin b.dat $size $nrhs
"$DIR_EXE"/$mat r.bin r.dat $size $nrhs
echo a
"$DIR_EXE"/$coo a.dat $size $size
echo b
"$DIR_EXE"/$coo b.dat $size $nrhs
echo lu
"$DIR_EXE"/$coo lu.dat $size $size
echo r
"$DIR_EXE"/$coo r.dat $size $nrhs

echo ""
echo ""
mpirun -n $np "$DIR_EXE"/inv $size
"$DIR_EXE"/$mat a.bin a.dat $size $size
"$DIR_EXE"/$mat lu.bin lu.dat $size $size
echo a
"$DIR_EXE"/$coo a.dat $size $size
echo lu
"$DIR_EXE"/$coo lu.dat $size $size
