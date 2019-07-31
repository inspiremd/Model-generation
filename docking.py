"""INSPIRE Docking - Put a ligand (SMILES) in a protein (PDB)

Usage:
  docking.py -s=<SMILES> -p=<PDB> -o=<PATH>
  docking.py (-h | --help)
  docking.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  -s=<SMILES>   Path to file containing SMILES representation of small molecule ligand.
  -p=<PDB>      Path to PDB file containing protein target for docking.
  -o=<PATH>     Path to write output.

"""
from docopt import docopt
from impress_md import interface_functions

if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE Docking 0.0.1')
    smiles = arguments['-s']
    pdb = arguments['-p']
    path = arguments['-o']
    interface_functions.RunDocking(smiles,pdb,path)
    interface_functions.ParameterizeSystem(path)



