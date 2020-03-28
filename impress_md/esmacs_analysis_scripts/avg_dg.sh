#/bin/bash
# $1 is the first command line argument which is the number of replicas.
if [ $# -ne 1 ];then
    echo "Incorrect no of arguments"
    exit
else
    num_reps=$1
fi

rm -f $PWD/*-rep-all.dat 
kbt=0.5961 # kcal/mol at 300 K
beta=1.6775
dir=$SCRIPTS_PATH
for ((i=1; i<=num_reps; i++)); do
check=$(for file in `ls $PWD/rep$i-*`; do if [ ! -s $file ]; then echo "$file empty"; fi; done | wc | awk '{print $1}')
if [ $check = 0 ]; then
    paste $PWD/rep$i-pb-com.dat $PWD/rep$i-pb-rec.dat $PWD/rep$i-pb-lig.dat | awk '{print $6-$12-$18}' | awk -f $dir/ave.awk >> $PWD/dg-tot-pb-rep-all.dat
    paste $PWD/rep$i-gb-com.dat $PWD/rep$i-gb-rec.dat $PWD/rep$i-gb-lig.dat | awk '{print $3+$2-$9-$8-$15-$14}' > $PWD/rep$i-Eint.dat
    Eint_mean=$(awk '{print $1}' $PWD/rep$i-Eint.dat | awk -f $dir/ave.awk | awk '{print $1}')
    awk -v Eint_mean=$Eint_mean -v beta=$beta -v kbt=$kbt 'BEGIN {val=0} {val += exp(beta*($1-Eint_mean))} END {print kbt*log(val/NR)}' $PWD/rep$i-Eint.dat >> $PWD/tds-var-rep-all.dat
fi
done

paste $PWD/dg-tot-pb-rep-all.dat $PWD/tds-var-rep-all.dat | awk '{print $1 + $3}' > $PWD/dg-tot-theo-pb-rep-all.dat
rm $PWD/rep*.dat $PWD/dg-tot-pb-rep-all.dat $PWD/tds-var-rep-all.dat
