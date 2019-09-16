"""INSPIRE Absolute BFE Calculator - Computed binding free energy estimates with alchemical methods

Usage:
  alchemy.py -i=<STRUCTURES> -n=<NANOSECONDS> -l=<NUM_LAMBDA>
  alchemy.py (-h | --help)
  alchemy.py --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  -i=<STRUCTURES>   Path to files containing AMBER format topology and coordinates of the protein, ligand and complex.
  -n=<NANOSECONDS>  Trajectory length of simulation to run.
  -l=<NUM_LAMBDA>   Number of interpolation points between the ligand-protein complex and the protein alone.
"""
# TODO: needs a way to input the number of lambdas/array and the ns of simulation
from docopt import docopt
from impress_md import interface_functions
import timeit
start = timeit.default_timer()

if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE Alchemical Calculator 0.0.1')
    path = arguments['-i']
    nsteps_per_iter = 1000 # 2 ps
    niter = round(500 * float(arguments['-n']))
    nlambda = arguments['-l']
    
    interface_functions.RunAlchemy(path,niter,nsteps_per_iter,nlambda)
    with open(f'{path}/alchemical.log',"w+") as logf:
        logf.write("Execution time (sec): {}\n".format(timeit.default_timer() - start))


