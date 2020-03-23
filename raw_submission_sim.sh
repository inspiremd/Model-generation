#!/bin/bash
 #BSUB -P CHM155_001
 #BSUB -W 2:00
 #BSUB -nnodes 4
 #BSUB -alloc_flags "gpudefault smt1"
 #BSUB -J Sim
 #BSUB -o Sim.o%J
 #BSUB -e Sim.e%J

 export OMP_NUM_THREADS=1
 COMP="com"
 INPATH="$MEMBERWORK/chm155/inpath/$COMP"
 source ~/.bash_profile
 source ~/.bashrc
 conda activate openmm
 module load cuda gcc spectrum-mpi
 cd $INPATH
 for i in {1..24};
 do
	 OUTPATH="$INPATH/rep$i"
	 mkdir -p $OUTPATH
	 cd $OUTPATH
	 rm -f traj.dcd sim.log
	 jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 python $INPATH/sim.py &
 done

 wait
