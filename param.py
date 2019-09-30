"""INSPIRE Parameterization - Parameterize a docked ligand-protein system

Usage:
  param.py -i=<PATH> [-a]
  param.py (-h | --help)
  param.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  -i=<PATH>     Path to directory containing a PDBs of the ligand, receptor, and complex.
  -a            Parameterize the docked system with AMBER [default uses OpenEye].

"""
from docopt import docopt
from impress_md import interface_functions
import timeit
start = timeit.default_timer()

if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE Param 0.0.1')
    path = arguments['-i']
    
    if arguments['-a']:
        interface_functions.ParameterizeAMBER(path)
    else:
        interface_functions.ParameterizeOE(path)

with open(f'{path}/param.log',"w+") as logf:
    logf.write("Param time (sec): {}\n".format(timeit.default_timer() - start))

