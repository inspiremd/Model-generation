# Model-generation
Python scripts to generate an MD-ready model from smiles strings and run simple free energy calculations 

There are two command line executables.

## 0_docking.py
* Takes a smiles and pdb, generates conformers, docks, scores, and parameterizes the ligand. Also writes the protein structure to a separate file (apo).
* The output is a set of simulation-ready structures (ligand, apo, and complex) and a file called metrics.csv, which has the docking score and associated uncertainties. Most uncertainties are 0 right now. There are other auxiliary files that are saved in the output directory.
* Dependencies: OpenEye, Ambertools, Ambermini, docopt

## 1_mmgbsa.py
* This command should only be run after docking.py. 
* Input is a path to the directory where the input coordinates and parameters are saved. This should be the output path from the docking.py command.
* Also takes the nanosecond length of the simulation. 0 corresponds to an energy minimization.
* Output adds to the metrics.csv file
* Dependencies: OpenMM, numpy, pymbar, docopt


To get the three metrics for a smiles, including a 5 ns simulation, pick a smiles and call
~~~
python 0_docking.py -s $SMILES -i "input/receptor.oeb" -o "test" -p
python 1_mmgbsa.py -p "test" -n 0
python 1_mmgbsa.py -p "test" -n 5
~~~
