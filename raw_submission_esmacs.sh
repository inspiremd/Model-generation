#!/bin/bash
 #BSUB -P CHM155_001
 #BSUB -W 0:10
 #BSUB -nnodes 24
 #BSUB -J MMPBSA
 #BSUB -o MMPBSA.o%J
 #BSUB -e MMPBSA.e%J
 #BSUB -w 977752

 export OMP_NUM_THREADS=1
 module load cuda/10.1.243 gcc/6.4.0 spectrum-mpi/10.3.1.2-20200121
 export CUDA_HOME=/sw/summit/cuda/10.1.243
 INPATH="$MEMBERWORK/chm155/inpath"
 source /gpfs/alpine/scratch/apbhati/chm155/AmberTools19/amber18/amber.sh 
 export LD_LIBRARY_PATH="/sw/summit/cuda/10.1.243/lib:${LD_LIBRARY_PATH}"
 cd $INPATH
 for i in {1..24};
 do
	 OUTPATH="$INPATH/rep$i"
	 cd $OUTPATH
         jsrun -n 1 -a 40 -c 40 -g 0 MMPBSA.py.MPI -O -i $INPATH/mmpbsa.in -sp $INPATH/complex.prmtop -cp $INPATH/com.prmtop -rp $INPATH/apo.prmtop -lp $INPATH/lig.prmtop -y traj.dcd > mmpbsa.log &
 done

 wait

 for i in {1..24}; do
    cd $INPATH/rep$i
    cat _MMPBSA_complex_gb.mdout.{0..39} > _MMPBSA_complex_gb.mdout.all
    cat _MMPBSA_complex_gb_surf.dat.{0..39} > _MMPBSA_complex_gb_surf.dat.all
    cat _MMPBSA_complex_pb.mdout.{0..39} > _MMPBSA_complex_pb.mdout.all
    cat _MMPBSA_ligand_gb.mdout.{0..39} > _MMPBSA_ligand_gb.mdout.all
    cat _MMPBSA_ligand_gb_surf.dat.{0..39} > _MMPBSA_ligand_gb_surf.dat.all
    cat _MMPBSA_ligand_pb.mdout.{0..39} > _MMPBSA_ligand_pb.mdout.all
    cat _MMPBSA_receptor_gb.mdout.{0..39} > _MMPBSA_receptor_gb.mdout.all
    cat _MMPBSA_receptor_gb_surf.dat.{0..39} > _MMPBSA_receptor_gb_surf.dat.all
    cat _MMPBSA_receptor_pb.mdout.{0..39} > _MMPBSA_receptor_pb.mdout.all
    rm _MMPBSA_*.{0..39} reference.frc *.pdb *.inpcrd *.mdin* *.out
 done

