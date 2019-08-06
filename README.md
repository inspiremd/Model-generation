# Model-generation
Python scripts to generate an MD-ready model from smiles strings

There are two command line executables.

## docking.py
* Takes a smiles and pdb, generates conformers, docks, scores, and parameterizes the ligand.
* The output is a set of simulation-ready structures (ligand, apo, and complex) and a file called metrics.csv, which has the docking score and associated uncertainties. Most uncertainties are 0 right now. There are other auxiliary files that are saved in the output directory.
* Dependencies: OpenEye, Ambertools, Ambermini

## mmgbsa.py
* This command should only be run after docking.py. 
* Takes the path to a directory where the input coordinates and parameters are saved. This should be the output path from the docking.py command.
* Also takes the nanosecond length of the simulation. 0 corresponds to an energy minimization
* Output adds to the metrics.csv file
* Dependencies: OpenMM, numpy, pymbar


To get all three metrics for a smiles, including a 1 ns mmgbsa run, pick a smiles and call
~~~
python docking.py -s SMILES -p "input/com_axitinib.pdb" -o "examples/example_2"
python mmgbsa.py -p "examples/example_2" -n 0
python mmgbsa.py -p "examples/example_2" -n 1 
~~~


Note
Antechamber is taking a long time (~10 minutes) to run on my computer here, as is OpenMM. I found it took long on tetrazole rings. Dave had issues with phosphates. 

## Todo
  - Make docking.py and mmgbsa.py time the simulation and output to a log file (also output any errors?)
  - In mmgbsa.py -- 'Non-optimal GB parameters detected for GB model %s' % gbmodel... Was this dealt with by using mbondi3 in tleap?

  - Make docking.py accept a receptor.oeb file ~ should be working
