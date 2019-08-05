#!/bin/bash
if [ $# -ne 8 ]
then
        echo "Wrong number of parameters"
        exit 1
fi

machine=${1}
app=${2}
exe=${3}
nbnode=${4}
nbcore=${5}
size=${6}
lang=${7}
rf=${8}

d1=$(date +%s.%N)
t=( $(mpirun -n $nbcore $exe -s $size) )
d2=$(date +%s.%N)
t_app=$(echo "$d2 $d1" | awk '{printf "%f", $1 - $2}')
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
{"machine":"$machine",\
"nb_cores":"$nbcore",\
"nb_nodes":"$nbnode",\
"test":"$app",\
"lang":"$lang",\
"datasize":"$size",\
"date":"$date",\
"time_io":"${t[2]}",\
"time_i":"${t[1]}",\
"time_calc":"${t[0]}",\
"time_app":"$t_app",\
"success":"$success",\
"comment":"Compilation -O2"}
EOF
} | tee -a $rf

