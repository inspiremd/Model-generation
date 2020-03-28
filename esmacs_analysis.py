"""INSPIRE ESMACS Analyser 

Usage:
  esmacs_analysis.py -i=<ENERGY_FILES> [-r=<REPLICAS>] [-t=<TYPE>] 
  esmacs_analysis.py (-h | --help)
  esmacs_analysis.py --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  -i=<ENERGY_FILES> Path to files containing output of esmacs.py. Should be the same as outpath for sim_esmacs.py 
  -r=<REPLICAS>     Ensemble size for ESMACAS [default=24].
  -t=<TYPE>         Type of ESMACS: 1-traj -> 'one', 2-traj -> 'two', 3-traj -> 'three' [default='one']. 

"""
import sys, os
import numpy as np
from docopt import docopt
from impress_md import interface_functions
import timeit
start = timeit.default_timer()

if __name__ == '__main__':
    arguments = docopt(__doc__, version='INSPIRE ESMACS Analyser 0.0.1')
    path = arguments['-i']
    if arguments['-t'] is None:
        typ = 'one'
    else:
        typ = arguments['-t']

    if arguments['-r'] is None:
        reps = 24
    else:
        reps = arguments['-r']

    scripts_path = os.path.join(str(os.getcwd()),'impress_md','esmacs_analysis_scripts')
    import subprocess
    with interface_functions.working_directory(path):
        subprocess.check_output(f'export SCRIPTS_PATH={scripts_path}', shell=True)
        subprocess.check_output(f'export PATH=$PATH:$SCRIPTS_PATH', shell=True)
    if typ == 'one':
        inpath = os.path.join(path,'com')
        with interface_functions.working_directory(inpath):
            subprocess.check_output(f'1traj.com {reps}', shell=True)
        newpath = os.path.join(inpath,'1-traj_analysis')
        with interface_functions.working_directory(newpath):
            subprocess.check_output(f'avg_dg.sh {reps}', shell=True)
            subprocess.check_output(f'bootstrap_1traj.sh {reps}', shell=True)
        dat = []
        final = open(f'{newpath}/resample-avg-theo-pb.dat', 'r')
        for line in final.readlines():
            dat.append(float(line.strip()))
        dg_avg = np.mean(dat)
        dg_err = np.std(dat)
        print("ESMACS (1-traj) free energy estimate is {}({})\n".format(dg_avg, dg_err))
        with open(f'{path}/metrics.csv','r') as metrics:
            dat = metrics.readlines()
        with open(f'{path}/metrics.csv','w') as metrics:
            metrics.write(dat[0].replace('\n',',ESMACS-1traj,ESMACS_err\n'))
            metrics.write(dat[1].replace('\n',',{},{}\n'.format(dg_avg,dg_err)))
    elif typ == 'two':
        inpath = os.path.join(path,'com')
        with interface_functions.working_directory(inpath):
            subprocess.check_output(f'2traj.com {reps}', shell=True)
        newpath = os.path.join(inpath,'2-traj_analysis')
        inpath = os.path.join(path,'lig')
        with interface_functions.working_directory(inpath):
            subprocess.check_output(f'mtraj_lig.com {reps}', shell=True)
            subprocess.check_output(f'mv m-traj_analysis/rep*-lig.dat {newpath}/', shell=True)
        with interface_functions.working_directory(newpath):
            subprocess.check_output(f'2traj_dg.sh {reps}', shell=True)
            subprocess.check_output(f'mtraj_g_lig.sh {reps}', shell=True)
            subprocess.check_output(f'bootstrap_2traj.sh {reps}', shell=True)
        dat = []
        final = open(f'{newpath}/resample-avg-pb.dat', 'r')  # temporarily without entropy
        for line in final.readlines():
            dat.append(float(line.strip()))
        dg_avg = np.mean(dat)
        dg_err = np.std(dat)
        print("ESMACS (2-traj) free energy estimate is {}({})\n".format(dg_avg, dg_err))
        with open(f'{path}/metrics.csv','r') as metrics:
            dat = metrics.readlines()
        with open(f'{path}/metrics.csv','w') as metrics:
            metrics.write(dat[0].replace('\n',',ESMACS-2traj,ESMACS_err\n'))
            metrics.write(dat[1].replace('\n',',{},{}\n'.format(dg_avg,dg_err)))
    elif typ == 'three':
        inpath = os.path.join(path,'com')
        with interface_functions.working_directory(inpath):
            subprocess.check_output(f'3traj.com {reps}', shell=True)
        newpath = os.path.join(inpath,'3-traj_analysis')
        inpath = os.path.join(path,'apo')
        with interface_functions.working_directory(inpath):
            subprocess.check_output(f'mtraj_rec.com {reps}', shell=True)
            subprocess.check_output(f'mv m-traj_analysis/rep*-rec.dat {newpath}/', shell=True)
        inpath = os.path.join(path,'lig')
        with interface_functions.working_directory(inpath):
            subprocess.check_output(f'mtraj_lig.com {reps}', shell=True)
            subprocess.check_output(f'mv m-traj_analysis/rep*-lig.dat {newpath}/', shell=True)
        with interface_functions.working_directory(newpath):
            subprocess.check_output(f'3traj_g_com.sh {reps}', shell=True)
            subprocess.check_output(f'mtraj_g_rec.sh {reps}', shell=True)
            subprocess.check_output(f'mtraj_g_lig.sh {reps}', shell=True)
            subprocess.check_output(f'bootstrap_3traj.sh {reps}', shell=True)
        dat = []
        final = open(f'{newpath}/resample-avg-pb.dat', 'r')    # temporarily without entropy
        for line in final.readlines():
            dat.append(float(line.strip()))
        dg_avg = np.mean(dat)
        dg_err = np.std(dat)
        print("ESMACS (3-traj) free energy estimate is {}({})\n".format(dg_avg, dg_err))
        with open(f'{path}/metrics.csv','r') as metrics:
            dat = metrics.readlines()
        with open(f'{path}/metrics.csv','w') as metrics:
            metrics.write(dat[0].replace('\n',',ESMACS-3traj,ESMACS_err\n'))
            metrics.write(dat[1].replace('\n',',{},{}\n'.format(dg_avg,dg_err)))
    else:
        sys.exit("Invalid value for argument TYPE (-t). Should be one, two or three only! Exiting.")
           
    with open(f'{path}/esmacs_analysis.log',"w+") as logf:
        logf.write("Execution time (sec): {}\n".format(timeit.default_timer() - start))



