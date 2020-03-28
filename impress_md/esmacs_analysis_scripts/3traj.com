#!/bin/bash
# $1 is the first command line argument which is the number of replicas.
if [ $# -ne 1 ];then
    echo "Incorrect no of arguments"
    exit
else
    num_reps=$1
fi

  dir_data=$PWD
  mkdir $dir_data/3-traj_analysis
  dir_out=$dir_data/3-traj_analysis
  for ((n=1; n<=num_reps; n++)); do
    dir_in=$dir_data/rep$n
    rm -f $dir_in/_MMPBSA_receptor_*.all
    rm -f $dir_in/_MMPBSA_ligand_*.all
    awk '($1=="BOND"||$1=="VDWAALS"||$1=="1-4")' $dir_in/_MMPBSA_complex_pb.mdout.all | cut -c11-24,36-49,64-77 | awk -f $SCRIPTS_PATH/t-ser-mmpb.awk > $dir_out/mmpb_com.dat
    awk '($1=="BOND"||$1=="VDWAALS"||$1=="1-4")' $dir_in/_MMPBSA_complex_gb.mdout.all | cut -c11-24,36-49,64-77 | awk -f $SCRIPTS_PATH/t-ser-mmpb.awk > $dir_out/mmgb_com.dat
    awk '($1 != "#Frame"){print $2*0.0072}' $dir_in/_MMPBSA_complex_gb_surf.dat.all > $dir_out/esurf_com.dat
    paste $dir_out/mmpb_com.dat $dir_out/esurf_com.dat | awk -f $SCRIPTS_PATH/t-ser-mmpbsa.awk > $dir_out/rep$n-pb-com.dat 
    paste $dir_out/mmgb_com.dat $dir_out/esurf_com.dat | awk -f $SCRIPTS_PATH/t-ser-mmgbsa.awk > $dir_out/rep$n-gb-com.dat 
  done

rm $dir_out/mm?b_???.dat $dir_out/esurf_???.dat
