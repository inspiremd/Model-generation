"""INSPIRE MMGBSA Calculator - Computed bind free energy estimates *fairly* cheaply.

Usage:
  mmgbsa.py -p=<STRUCTURES> [-n=<NANOSECONDS>] [-r=<REPLICAS>] [-m]
  mmgbsa.py (-h | --help)
  mmgbsa.py --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  -p=<STRUCTURES>   Path to files containing AMBER format topology and coordinates of the protein, ligand and complex.
  -n=<NANOSECONDS>  Trajectory length of simulation to run (0 = minimization alone) [default=0].
  -r=<REPLICAS>     Number of replica simulations to execute [default=1]. NOT CURRENTLY IN USE
  -m                Minimize components independently. CURRENTLY ALWAYS TRUE

"""
from docopt import docopt
from impress_md import interface_functions
import timeit
start = timeit.default_timer()

if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE MMGBSA Calculator 0.0.1')
    path = arguments['-p']
    one_traj = arguments['-m']

    if arguments['-n'] is None or float(arguments['-n'])*1000 == 0:
        interface_functions.RunMinimization(path, path, one_traj)
        with open(f'{path}/minimization.log',"w+") as logf:
            logf.write("Execution time (sec): {}\n".format(timeit.default_timer() - start))
    else:
        niter = round(float(arguments['-n'])*1000)
        interface_functions.RunMMGBSA(path, path, niter)
        with open(f'{path}/mmgbsa.log',"w+") as logf:
            logf.write("Execution time (sec): {}\n".format(timeit.default_timer() - start))



