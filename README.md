# Model-generation
Python scripts to generate an MD-ready model from smiles strings

There are three command line executables.

## 0_docking.py
* Takes a smiles and pdb, generates conformers, docks, scores, and parameterizes the ligand.
* The output is a set of simulation-ready structures (ligand, apo, and complex) and a file called metrics.csv, which has the docking score and associated uncertainties. Most uncertainties are 0 right now. There are other auxiliary files that are saved in the output directory.
* Dependencies: OpenEye

## 1_ante.sh
* Input a directory where the results of the docking step was saved
* Antechamber is run to prepare the ligand for simulation
* Dependencies: Ambertools, Ambermini

## 2_mmgbsa.py
* This command should only be run after docking.py. 
* Takes the path to a directory where the input coordinates and parameters are saved. This should be the output path from the docking.py command.
* Also takes the nanosecond length of the simulation. 0 corresponds to an energy minimization. <b> This does not work properly for n != 0.</b>
* Output adds to the metrics.csv file
* Dependencies: OpenMM, numpy, pymbar


To get the first two metrics for a smiles, pick a smiles and call
~~~
python 0_docking.py -s $SMILES -i "input/com_axitinib.pdb" -o "test"
./1_ante.sh test
python 2_mmgbsa.py -p "test" -n 0 
~~~


Note
Antechamber is taking a long time (~10 minutes) to run on my computer here, as is OpenMM. I found it took long on tetrazole rings. Dave had issues with phosphates & things that should be charged. 
Something is going wrong with 2_mmgbsa for -n 1. It runs locally and on Summit's login nodes, but when submitted through jsrun.
