"""INSPIRE Docking - Put a ligand (SMILES) in a protein (PDB)

Usage:
  docking.py -s=<SMILES> -o=<PATH> -i=<PATH>
  docking.py (-h | --help)
  docking.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  -s=<SMILES>   Path to file containing SMILES representation of small molecule ligand.
  -i=<PATH>     Path to PDB or OEB file containing protein target for docking.
                The openeye binary (OEB) can be precompiled to speed up docking time.
  -o=<PATH>     Path to write output.
"""
from docopt import docopt
from impress_md import interface_functions
import timeit
start = timeit.default_timer()

if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE Docking 0.0.1')
    smiles = arguments['-s']
    struct = arguments['-i']
    
    path = arguments['-o']

    interface_functions.RunDocking(smiles,struct,path)
    docked_time = timeit.default_timer()

with open(f'{path}/docking.log',"w+") as logf:
    logf.write("Docking time (sec): {}\n".format(docked_time - start))

