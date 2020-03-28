"""INSPIRE ESMACS Calculator - Computes binding free energy estimates taking conformations from MD simulation(s).

Usage:
  esmacs.py -i=<STRUCTURES> [-o=<TRAJECTORY>] [-r=<REPLICAS>] [-c=<COMPONENT>]
  esmacs.py (-h | --help)
  esmacs.py --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  -i=<STRUCTURES>   Path to input files containing AMBER format topology of the system (same as the input path used for sim_esmacs.py).
  -o=<TRAJECTORY>   Path to files containing output of sim_esmacs.py [default=same as input path].
  -r=<REPLICAS>     Ensemble size for ESMACS; same as that used in sim_esmacs.py NOT CURRENTLY IN USE (To be worked out with RCT team).
  -c=<COMPONENT>    complex -> 'com', protein -> 'apo', ligand -> 'lig' [default='com']. 'com' is always needed; others only in case of 2- or 3-traj ESMACS when independent simulations have been run for other components.

"""
import sys, os
from docopt import docopt
from impress_md import interface_functions
import timeit
start = timeit.default_timer()

if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE ESMACS Calculator 0.0.1')
    inpath = arguments['-i']
#    replicas = arguments['-r'] # to be worked out with RCT team
    if arguments['-o'] is None:
        outpath = str(inpath)
    else:
        outpath = arguments['-o']

    if arguments['-c'] is None:
        comp = 'com'
    else:
        comp = arguments['-c']

    path = os.getcwd
    outpath1 = os.path.join(outpath,comp,'rep1') # rep1 is temporary; to be worked out with RCT team
    import subprocess
    with interface_functions.working_directory('outpath1'):
        if comp == 'com':
            subprocess.check_output(f'MMPBSA.py.MPI -O -i {path}/impress_md/esmacs.in -sp {inpath}/{comp}_sol.prmtop -cp {inpath}/com.prmtop -rp {inpath}/apo.prmtop -lp {inpath}/lig.prmtop -y {outpath1}/traj.dcd > {outpath1}/mmpbsa.log',shell=True) # number of MPI ranks = number of ns (excluding first 2 ns) * 10; 40 is default
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_complex_gb.mdout.* | sort -V | xargs cat > {outpath}/_MMPBSA_complex_gb.mdout.all', shell=True)
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_complex_gb_surf.dat.* | sort -V | xargs cat > {outpath}/_MMPBSA_complex_gb_surf.dat.all', shell=True)
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_complex_pb.mdout.* | sort -V | xargs cat > {outpath}/_MMPBSA_complex_pb.mdout.all', shell=True)
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_ligand_gb.mdout.* | sort -V | xargs cat > {outpath}/_MMPBSA_ligand_gb.mdout.all', shell=True)
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_ligand_gb_surf.dat.* | sort -V | xargs cat > {outpath}/_MMPBSA_ligand_gb_surf.dat.all', shell=True)
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_ligand_pb.mdout.* | sort -V | xargs cat > {outpath}/_MMPBSA_ligand_pb.mdout.all', shell=True)
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_receptor_gb.mdout.* | sort -V | xargs cat > {outpath}/_MMPBSA_receptor_gb.mdout.all', shell=True)
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_receptor_gb_surf.dat.* | sort -V | xargs cat > {outpath}/_MMPBSA_receptor_gb_surf.dat.all', shell=True)
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_receptor_pb.mdout.* | sort -V | xargs cat > {outpath}/_MMPBSA_receptor_pb.mdout.all', shell=True)
            subprocess.check_output(f'rm {outpath1}/_MMPBSA_*.[0-9]* {outpath1}/reference.frc {outpath1}/*.inpcrd {outpath1}/*.mdin* {outpath1}/*.out', shell=True)
         else:
            subprocess.check_output(f'MMPBSA.py.MPI -O -i {path}/impress_md/esmacs.in -sp {inpath}/{comp}_sol.prmtop -cp {inpath}/{comp}.prmtop -y {outpath1}/traj.dcd > {outpath1}/mmpbsa.log',shell=True) # number of MPI ranks = number of ns (excluding first 2 ns) * 10; 40 is default
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_complex_gb.mdout.* | sort -V | xargs cat > {outpath}/_MMPBSA_complex_gb.mdout.all', shell=True)
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_complex_gb_surf.dat.* | sort -V | xargs cat > {outpath}/_MMPBSA_complex_gb_surf.dat.all', shell=True)
            subprocess.check_output(f'ls {outpath1}/_MMPBSA_complex_pb.mdout.* | sort -V | xargs cat > {outpath}/_MMPBSA_complex_pb.mdout.all', shell=True)
            subprocess.check_output(f'rm {outpath1}/_MMPBSA_*.[0-9]* {outpath1}/reference.frc {outpath1}/*.inpcrd {outpath1}/*.mdin* {outpath1}/*.out', shell=True)

    with open(f'{outpath}/esmacs.log',"w+") as logf:
        logf.write("Execution time (sec): {}\n".format(timeit.default_timer() - start))



