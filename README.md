This runs of x86 ONLY for the moment. 

Build the conda enviroment using the yml file or just run stuff and figure out how to make it work (it's not rocket science ;)

You can run end to end with 
```python 
mpiexec -np 8 python runner.py input/john_smiles_kinasei.smi
```

You can alter how those runs peform using the policy file.

In the future, this workflow will be moved to Radical/ENTK, but for now this works. 



---- 

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

## sim.py
* Performs MD simulation of the given system after solvating it with explicit water molecules.
* This command should only be run after docking.py and param.py.
* Input is path to the directory where the coordinates and parameters (AMBER format) of the system are saved. This should be the output path for docking.py.
* Takes as input the component to be simulated (complex, protein or ligand) and the nanosecond length of the simulation. 0 corresponds to energy minimisation.
* Output is added to the metrics.csv file.
* Dependencies: OpenMM, docopt

## sim_esmacs.py
* Performs MD simulation of the given solvated system and saves trajectories for ESMACS calculations.
* Should only be run after docking and parameterisation steps have been successful.
* Input path is path to the directory containing the AMBER topology and coordinate files (output path for `docking.py`).
* Other inputs: component to be simulated (com, apo or lig), nanosecond length of simulation (first 2 ns considered equilibration) and ensemble size for ESMACS.
* Should be called separately for each component in case of multiple trajectory version of ESMACS.
* Dependencies: OpenMM, docopt

## mmgbsa_explicit.py
* Computes binding free energy estimates (MMGBSA) taking conformations from MD simulations(s).
* This command should only be run after `sim.py`.
* Input path is the location of parameters (AMBER format) for the system. Should be the same as used for `sim.py`.
* Output path is the location of trajectory files, that is the output path of `sim.py`.
* User has the choice of performing 1-traj or 3-traj versions of MMGBSA.
* Output is added to the `metrics.csv` file.
* Dependencies: OpenMM, docopt

## esmacs.py
* Performs ESMACS related calculations on the solvated trajectories.
* Should be only be executed after `sim_esmacs.py` has been successful.
* Input and output paths should be the same as those used for `sim_esmacs.py`.
* Other inputs: Ensemble size for ESMACS, component to be analysed.
* Should be called separately for each component in case of multiple trajectory version of ESMACS.
* Dependencies: AmberTools19, docopt

## esmacs_analysi.py
* Performs final ESMACS analysis on the MMPBSA output of AmberTools to give binding affinity along with its statistial uncertainty.
* Should only be run after successfuly execution of `esamcs.py`.
* Can be run on login nodes of Summit/local machine.
* Input path should be the output path of `sim_esmacs.py` and `esmacs.py`.
* Other inputs: Ensemble size for ESMACS, type of calculation (1, 2 or 3 traj).
* Dependencies: docopt

## alchem.py
* Uses an alchemical method to calculate the absolute binding free energy of a ligand.
* User enters the number of lambda windows and the length of simulation at each window.
* The windows run in series but this could be parallelized.
* Applying constraints should increase the convergence of the system

To get the four metrics for a smiles, including a 5 ns simulation, pick a smiles and call
~~~bash
python docking.py -s $SMILES -i "input/receptor.oeb" -o "test"
python param.py -i "test"
python mmgbsa.py -p "test" -n 0
python mmgbsa.py -p "test" -n 5
python alchem.py -i "test" -l 6 -n 5
~~~
