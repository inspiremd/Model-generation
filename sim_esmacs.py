"""INSPIRE ESMACS Simulator 

Usage:
  sim_esmacs.py -i=<STRUCTURES> [-o=<TRAJECTORY>] [-c=<COMPONENT>] [-n=<NANOSECONDS>] [-r=<REPLICAS>]
  sim_esmacs.py (-h | --help)
  sim_esmacs.py --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  -i=<STRUCTURES>   Path to input files containing AMBER format topology and coordinates of the solvated system.
  -o=<TRAJECTORY>   Path to output files (trajectory, energy, temp, etc.) of the simulation [default=same as input path].
  -c=<COMPONENT>    Component (complex=com, protein=apo and ligand=lig) which needs to be simulated [default=com]. For 1-traj ESMACS, only 'com' needs to be simulated but for multiple traj versions, other components need to be simulated too.
  -n=<NANOSECONDS>  Length of simulation to run in nanoseconds [default=6]. First 2 ns are ignored as equilibration and conformations are taken only from the remainder. So, it should always be greater than 2.
  -r=<REPLICAS>     Ensemble size for ESMACS (Standard is 25. Better be a multiple of 6 on Summit; so 24 preferable). Ignored here; Implementation to be worked out with RCT team.
"""
from docopt import docopt
from impress_md import interface_functions
import timeit
start = timeit.default_timer()

if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE ESMACS Simulator 0.0.1')
    inpath = arguments['-i']
    comp = arguments['-c']
#    replicas = arguments['-r']
    if arguments['-o'] is None:
        outpath = str(inpath)
    else:
        outpath = arguments['-o']

    if arguments['-n'] is None:
        nsteps = 3000000
    else:
        nsteps = round(float(arguments['-n'])*500000)  # assuming timestep of 2 fs

    interface_functions.Simulation_ESMACS(inpath, outpath, nsteps, comp)
    with open(f'{outpath}/{comp}_simulation_ESMACS.log',"w+") as logf:
        logf.write("Execution time (sec): {}\n".format(timeit.default_timer() - start))


