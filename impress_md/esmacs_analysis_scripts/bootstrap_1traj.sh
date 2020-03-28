#/bin/bash
# $1 is the first command line argument which is the number of replicas.
if [ $# -ne 1 ];then
    echo "Incorrect no of arguments"
    exit
else
    num_reps=$1
fi

module load gcc/4.8.5
module load r/3.5.0
for i in theo-pb; do
	Rscript $SCRIPTS_PATH/bootstrap_1traj.R $i $num_reps
	rm $PWD/dg-tot-$i-rep-all.dat
done

