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

t=( $(mpirun -n $nbcore $exe -s $size) )
echo ${t[*]}
success="true"

if [ ${#t[*]} -ne 3 ]
then
	t="error"
	success="false"
fi

date=$(date +%Y%m%d-%H%M%S)

{
#echo "$machine;$nbhosts;$nbnodes;$app;$blocks;$size;$procs;$date;$nbWorker;$tmpGen;$totGen;$tmpNoG;$totNoG;$success;$nbTask"
cat << EOF
{"machine":"Poincare",\
"nb_cores":"$nbcore",\
"nb_nodes":"$nbnode",\
"test":"$app",\
"lang":"$lang",\
"datasize":"$size",\
"date":"$date",\
"time_io":"${t[2]}",\
"time_i":"${t[1]}",\
"time_calc":"${t[0]}",\
"success":"$success",\
"comment":"Compilation -O2"}
EOF
} | tee -a ~/results.csv

