"""INSPIRE MMGBSA Calculator - Computed bind free energy estimates *fairly* cheaply.

Usage:
  mmgbsa.py -t=<TOPOLOGY> -c=<COORDINATES> [-n=<NANOSECONDS>] [-r=<REPLICAS>] [-m]
  mmgbsa.py (-h | --help)
  mmgbsa.py --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  -t=<TOPOLOGY>     Path to file containing AMBER format topology of protein-ligand complex.
  -c=<COORDINATES>  Path to file containing AMBER format coordinates of protein-ligand complex.
  -n=<NANOSECONDS>  Trajectory length of simulation to run (0 = minimization alone) [default=0].
  -r=<REPLICAS>     Number of replica simulations to execute [default=1].
  -m                Minimize components independently.

"""
from docopt import docopt


if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE MMGBSA Calculator 0.0.1')
    print(arguments)