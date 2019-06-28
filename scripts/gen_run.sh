if [ $# -ne 3 ]
then
	echo "common use"
	echo "./compileAll.sh matSize processes nodes"
	echo "./compileAll.sh 16384 1024 64"
	echo "./compileAll.sh 160 16 1"
	exit
fi

sed -E "s/PROCS/$2/g;s/NODES/$3/g;s/SIZE/$1/g;" .submit_sca.sh > submit_sca_$1.sh
sed -E "s/PROCS/$2/g;s/NODES/$3/g;s/SIZE/$1/g;" .submit_mpi.sh > submit_mpi_$1.sh
sed -E "s/PROCS/$2/g;s/NODES/$3/g;s/SIZE/$1/g;" .submit_xmp.sh > submit_xmp_$1.sh
