#!/bin/bash
if [ $# -ne 6 ]
then
	echo "Wrong number of parameters"
	exit 1
fi

app=$1
size=$2
exe=$3
nbnode=$4
nbcore=$5
lang=$6

t=( $(mpirun -n $nbcore $exe $size) )
echo ${t[*]}
success="true"

if [ ${#t[*]} -ne 3 ]
then
	t="error"
	success="false"
fi

echo "Poincare;$nbcore;$nbnode;$app;${lang};$size;$nbcore;$(date);-1;0;0;${t[0]};0;$success;;Compilation -O2" >> ~/results.csv
echo "Poincare;$nbcore;$nbnode;$app;${lang}L;$size;$nbcore;$(date);-1;0;0;${t[1]};0;$success;;Compilation -O2" >> ~/results.csv
echo "Poincare;$nbcore;$nbnode;$app;${lang}LS;$size;$nbcore;$(date);-1;0;0;${t[2]};0;$success;;Compilation -O2" >> ~/results.csv

