#/bin/bash
# $1 is the first command line argument which is the number of replicas.
if [ $# -ne 1 ];then
    echo "Incorrect no of arguments"
    exit
else
    num_reps=$1
fi

rm -f *com*rep-all.dat 
dir=$SCRIPTS_PATH
for ((i=1; i<=num_reps; i++)); do
check=$(for file in `ls rep$i-*-com.dat`; do if [ ! -s $file ]; then echo "$file empty"; fi; done | wc | awk '{print $1}')
if [ $check = 0 ]; then
    awk '{print $6}' $PWD/rep$i-pb-com.dat | awk -f $dir/ave.awk >> $PWD/g-com-tot-pb-rep-all.dat
fi
done

rm $PWD/rep*com.dat
