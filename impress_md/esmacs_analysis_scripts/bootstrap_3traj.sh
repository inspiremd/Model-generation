#/bin/bash
# $1 is the first command line argument which is the number of replicas.
if [ $# -ne 1 ];then
    echo "Incorrect no of arguments"
    exit
else
    num_reps=$1
fi

rm -f dg-final.dat
module load gcc/4.8.5
module load r/3.5.0
for i in pb; do
	Rscript $SCRIPTS_PATH/bootstrap_3traj.R $i $num_reps
	rm $PWD/g-*-tot-$i-rep-all.dat
done

