"""INSPIRE MMGBSA Calculator - Computes binding free energy estimates taking conformations from MD simulation(s).

Usage:
  mmgbsa_explicit.py -i=<STRUCTURES> [-o=<TRAJECTORY>] [-r=<REPLICAS>] [-m]
  mmgbsa_explicit.py (-h | --help)
  mmgbsa_explicit.py --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  -i=<STRUCTURES>   Path to input files containing AMBER format topology of the system (same as the input path used for sim.py).
  -o=<TRAJECTORY>   Path to files containing output of sim.py (output of MD simulation(s) using explicit solvent) [default=same as input path].
  -r=<REPLICAS>     Number of replica simulations to execute [default=1]. NOT CURRENTLY IN USE
  -m                Use independent trajectory for each component. yes --> 3-traj; no --> 1-traj [default=no].

"""
import sys, os
from docopt import docopt
from impress_md import interface_functions
import timeit
start = timeit.default_timer()

if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE MMGBSA (explicit solvent) Calculator 0.0.1')
    inpath = arguments['-i']

    if arguments['-m'] is None:
        three_traj = 'no'
    else:
        three_traj = arguments['-m']

    if arguments['-o'] is None:
        outpath = str(inpath)
    else:
        outpath = arguments['-o']

    del_G = interface_functions.RunMMGBSA_explicit(inpath, outpath, three_traj)
    with open(f'{outpath}/mmgbsa_explicit.log',"w+") as logf:
        logf.write("The binding affinity (MMGBSA) of the system is {} kJ/mol\n".format(del_G))
        logf.write("Execution time (sec): {}\n".format(timeit.default_timer() - start))



