import sys, os
from contextlib import contextmanager

@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)

def RunDocking(smiles, inpath, outpath, padding=4):
    from . import conf_gen
    from . import dock_conf
    # if not os.path.exists(outpath):
    #     os.mkdir(outpath)
    confs = conf_gen.SelectEnantiomer(conf_gen.FromString(smiles))
    # This receptor can be pre-compiled to an oeb. It may speed things up notably
    filename, file_extension = os.path.splitext(inpath)
    if file_extension == ".oeb":
        receptor = dock_conf.PrepareReceptorFromBinary(inpath)
    else: # else it is a pdb
        receptor = dock_conf.PrepareReceptor(inpath,padding,outpath)
    dock, lig = dock_conf.DockConf(receptor,confs,MAX_POSES=1)

    dock_conf.WriteStructures(receptor, lig, f'{outpath}/apo.pdb', f'{outpath}/lig.pdb')
    with open(f'{outpath}/metrics.csv','w+') as metrics:
        metrics.write("Dock,Dock_U\n")
        metrics.write("{},{}\n".format(dock_conf.BestDockScore(dock,lig),0))
    from openeye import oedepict
    oedepict.OEPrepareDepiction(lig)
    oedepict.OERenderMolecule(f'{outpath}/lig.png',lig)

## Somehow hop over to the path
def ParameterizeSystem(path):
    import subprocess
    with working_directory(path):
        subprocess.check_output(f'antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y',shell=True)
        subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod',shell=True)
        # Wrap tleap to get $path/(com|lig|apo).(inpcrd|prmtop)
        with open(f'leap.in','w+') as leap:
            leap.write("source leaprc.protein.ff14SBonlysc\n")
            leap.write("source leaprc.gaff\n")
            leap.write("set default PBRadii mbondi3\n")
            leap.write("rec = loadPDB apo.pdb # May need full filepath?\n")
            leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")
            leap.write("lig = loadmol2 lig.mol2\n")
            leap.write("loadAmberParams lig.frcmod\n")
            leap.write("com = combine {rec lig}\n")
            leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
            leap.write("saveAmberParm com com.prmtop com.inpcrd\n")
            leap.write("quit\n")
        subprocess.check_output(f'tleap -f leap.in',shell=True)
    

def RunMinimization(build_path, outpath, one_traj=False):
    """
    We are minimizing all three structures, then checking the potential energy using GB forcefields
    We could, alternatively, minimize the docked structure and then extract trajectories (1 frame long), more like a 1-trajectory mmgbsa.
    output is path used in "RunDocking". It has a metric.csv file.
    Notes
        MinimizedEnergy function requires path/prefix to inpcrd and prmtop files.
        These are created by the ParameterizeSystem function...
    """
    from . import minimize
    success = True
    try:
        rec_energy = minimize.MinimizedEnergy(f'{build_path}/apo')
        lig_energy = minimize.MinimizedEnergy(f'{build_path}/lig')
        com_energy = minimize.MinimizedEnergy(f'{build_path}/com')
        diff_energy = com_energy - lig_energy - rec_energy
    except:
        success = False
    
    if one_traj:
        print("1-traj calculation not ready")

    with open(f'{outpath}/metrics.csv','r') as metrics:
        dat = metrics.readlines()
    with open(f'{outpath}/metrics.csv','w') as metrics:
        metrics.write(dat[0].replace('\n',',Minimize,Minimize_U\n'))
        if success:
            metrics.write(dat[1].replace('\n',',{},{}\n'.format(diff_energy,0)))
        else:
            metrics.write(dat[1].replace('\n',',NA,NA\n'))


def RunMMGBSA(inpath, outpath, niter=1000):
    from . import mmgbsa
    crds = {'lig':f'{inpath}/lig.inpcrd','apo':f'{inpath}/apo.inpcrd','com':f'{inpath}/com.inpcrd'}
    prms = {'lig':f'{inpath}/lig.prmtop','apo':f'{inpath}/apo.prmtop','com':f'{inpath}/com.prmtop'}
    enthalpies = mmgbsa.simulate(crds, prms, niter)
    mmgbsa.subsample(enthalpies)
    energies = mmgbsa.mmgbsa(enthalpies)
    # Now write this to file
    with open(f'{outpath}/metrics.csv','r') as metrics:
        dat = metrics.readlines()
    with open(f'{outpath}/metrics.csv','w') as metrics:
        metrics.write(dat[0].replace('\n',',mmgbsa,mmgbsa_U\n'))
        metrics.write(dat[1].replace('\n',',{},{}\n'.format(energies[0]['diff'],energies[1]['diff'])))
    return energies

def GetMetrics(smiles, pdb, outpath, method):
    """Method can be
        'dock'
        'minimize'
        'mmgbsa'
        'absolute' - Not in use
    """
    if method == 'dock':
        RunDocking(smiles,pdb,outpath)
        ParameterizeSystem(outpath)
        RunMinimization(outpath,outpath)
    if method == 'mmgbsa':
        RunMMGBSA(outpath, 500)
