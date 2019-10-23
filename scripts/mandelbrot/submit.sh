#!/bin/bash

if [ $# -ne 4 ]
then
        echo "Wrong number of parameters"
        exit 1
fi

timewall="00:30:00"
nodes=$1
cores=$(($nodes * 20))
datasize=$2
lang=$3
app=$4

lang_lower=$(echo ${lang} | tr '[:upper:]' '[:lower:]')
if [ "${app}" == "blockLU" ]; then exe="${lang_lower}_lu";
elif [ "${app}" == "blockGauss" ]; then exe="${lang_lower}_sls_g";
elif [ "${app}" == "blockGaussJordan" ]; then exe="${lang_lower}_sls_gj";
elif [ "${app}" == "blockLUsolveLS" ]; then exe="${lang_lower}_sls_lu";
elif [ "${app}" == "gaussJordan_inv" ]; then exe="${lang_lower}_inv_gj";
else echo "Wrong app name"; exit 1; fi

sed "s/TIMEWALL/${timewall}/;
     s/NNODES/${nodes}/;
     s/DATASIZE/${datasize}/g;
     s/APP/${app}/g;
     s%EXE%executables/${exe}%g;
     s/LANG/${lang}/g;
     s/NCORES/${cores}/" scripts/mandelbrot/launch_template.sh > submit_${lang}_${app}_${datasize}_${nodes}.sh

sh submit_${lang}_${app}_${datasize}_${nodes}.sh
