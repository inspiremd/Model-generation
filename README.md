# Model-generation
Python scripts to generate an MD-ready model from smiles strings

There are two command line executables.

## 0_docking.py
* Takes a smiles and pdb, generates conformers, docks, scores, and parameterizes the ligand.
* The output is a set of simulation-ready structures (ligand, apo, and complex) and a file called metrics.csv, which has the docking score and associated uncertainties. Most uncertainties are 0 right now. There are other auxiliary files that are saved in the output directory.
* Dependencies: OpenEye, Ambertools, Ambermini

## 1_mmgbsa.py
* This command should only be run after docking.py. 
* Takes the path to a directory where the input coordinates and parameters are saved. This should be the output path from the docking.py command.
* Also takes the nanosecond length of the simulation. 0 corresponds to an energy minimization. <b> This does not work properly for n != 0.</b>
* Output adds to the metrics.csv file
* Dependencies: OpenMM, numpy, pymbar


To get the first two metrics for a smiles, pick a smiles and call
~~~
python 0_docking.py -s $SMILES -i "input/receptor.oeb" -o "test" -p
python 1_mmgbsa.py -p "test" -n 0 
~~~


Note

Parameterization is done using OpenEye's AM1BCC charges. 
Something is going wrong with mmgbsa for -n 1. It runs locally and on Summit's login nodes, but CUDA is not properly recognized when submitted through jsrun.
