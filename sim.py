"""INSPIRE MD Simulator - Compute potential energy of the system solvated with explicit water.

Usage:
  sim.py -i=<STRUCTURES> [-o=<TRAJECTORY>] [-c=<COMPONENT>] [-n=<NANOSECONDS>]
  sim.py (-h | --help)
  sim.py --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  -i=<STRUCTURES>   Path to input files containing AMBER format topology and coordinates of the system.
  -o=<TRAJECTORY>   Path to output files (trajectory, energy, temp, etc.) of the simulation [default=same as input path].
  -c=<COMPONENT>    Component (complex=com, protein=apo and ligand=lig) which needs to be simulated [default=com].
  -n=<NANOSECONDS>  Length of simulation to run in nanoseconds (0 = minimization alone) [default=0].

"""
from docopt import docopt
from impress_md import interface_functions
import timeit
start = timeit.default_timer()

if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE MD Simulator 0.0.1')
    inpath = arguments['-i']
    comp = arguments['-c']
    if arguments['-o'] is None:
        outpath = str(inpath)
    else:
        outpath = arguments['-o']

    if arguments['-n'] is None:
        nsteps = 0
    else:
        nsteps = round(float(arguments['-n'])*500000)  # assuming timestep of 2 fs

    potential = interface_functions.Simulation_explicit(inpath, outpath, nsteps, comp)
    with open(f'{outpath}/simulation_explicit.log',"w+") as logf:
        logf.write("Potential energy of the simulated system is {} kJ/mol.\n".format(potential))
        logf.write("Execution time (sec): {}\n".format(timeit.default_timer() - start))


