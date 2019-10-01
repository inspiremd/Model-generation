# Model-generation
Python scripts to generate an MD-ready model from smiles strings and run simple free energy calculations 

There are two command line executables.

## docking.py
* Takes a smiles and pdb, generates conformers, docks, and scores the ligand.
* The output is a set of simulation-ready structures (ligand, apo, and complex) and a file called metrics.csv, which has the docking score and associated uncertainties. Most uncertainties are 0 right now. There are other auxiliary files that are saved in the output directory.
* Dependencies: OpenEye, Ambertools, Ambermini, docopt

## param.py
* Parameterizes the ligand using either OpenEye or Amber.

## mmgbsa.py
* This command should only be run after docking.py. 
* Input is a path to the directory where the input coordinates and parameters are saved. This should be the output path from the docking.py command.
* Also takes the nanosecond length of the simulation. 0 corresponds to an energy minimization.
* Output adds to the metrics.csv file
* Dependencies: OpenMM, numpy, pymbar, docopt

## alchem.py
* Uses an alchemical method to calculate the absolute binding free energy of a ligand.
* User enters

To get the four metrics for a smiles, including a 5 ns simulation, pick a smiles and call
~~~bash
python docking.py -s $SMILES -i "input/receptor.oeb" -o "test"
python param.py -i "test"
python mmgbsa.py -p "test" -n 0
python mmgbsa.py -p "test" -n 5
python alchem.py -i "test" -l 6 -n 5
~~~
