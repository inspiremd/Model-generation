import sys, os
import conf_gen, dock_conf


def RunDocking(smiles, pdb, outpath, padding=4):
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    confs = conf_gen.SelectEnantiomer(conf_gen.FromString(smiles))
    # This receptor can be pre-compiled to an oeb. It may speed things up notably
    if not os.path.exists(f'{outpath}/receptor.oeb'):
        receptor = dock_conf.PrepareReceptor(pdb,padding,outpath)
    else:
        receptor = dock_conf.PrepareReceptorFromBinary(f'{outpath}/receptor.oeb')
    dock, lig = dock_conf.DockConf(receptor,confs,MAX_POSES=1)

    dock_conf.WriteStructures(receptor, lig, f'{outpath}/apo.pdb', f'{outpath}/lig.pdb')
    with open(f'{outpath}/metrics.csv','w+') as metrics:
        metrics.write("Dock,Dock_U\n")
        metrics.write("{},{}\n".format(dock_conf.BestDockScore(dock,lig),0))


def ParameterizeSystem(path):
    import subprocess
    # Run antechamber at f'{path}/lig.pdb'
    subprocess.check_output(f'antechamber -i {path}/lig.pdb -fi pdb -o {path}/lig.mol2 -fo mol2 -c bcc -pf y -an y',shell=True)
    subprocess.check_output(f'parmchk2 -i {path}/lig.mol2 -f mol2 -o {path}/lig.frcmod',shell=True)
    subprocess.check_output('mv sqm* {}/'.format(path),shell=True)
    
    # Wrap tleap to get $path/(com|lig|apo).(inpcrd|prmtop)
    with open(f'{path}/leap.in','w+') as leap:
        leap.write("source leaprc.protein.ff14SBonlysc\n")
        leap.write("source leaprc.gaff\n")
        leap.write("set default PBRadii mbondi3\n")
        leap.write("rec = loadPDB {}/apo.pdb # May need full filepath?\n".format(path))
        leap.write("saveAmberParm rec {}/apo.prmtop {}/apo.inpcrd\n".format(path,path))
        leap.write("lig = loadmol2 {}/lig.mol2\n".format(path))
        leap.write("loadAmberParams {}/lig.frcmod\n".format(path))
        leap.write("com = combine {rec lig}\n")
        leap.write("saveAmberParm lig {}/lig.prmtop {}/lig.inpcrd\n".format(path,path))
        leap.write("saveAmberParm com {}/com.prmtop {}/com.inpcrd\n".format(path,path))
        leap.write("quit\n")
    subprocess.check_output(f'tleap -f {path}/leap.in',shell=True)
    subprocess.check_output('mv leap.log {}/'.format(path),shell=True)


def RunMinimization(build_path, outpath, three_traj=True):
    """
    We are minimizing all three structures, then checking the potential energy using GB forcefields
    We could, alternatively, minimize the docked structure and then extract trajectories (1 frame long), more like a 1-trajectory mmgbsa.
    output is path used in "RunDocking". It has a metric.csv file.
    Notes
        MinimizedEnergy function requires path/prefix to inpcrd and prmtop files.
        These are created by the ParameterizeSystem function...
    """
    from minimize import MinimizedEnergy
    success = True
    try:
        rec_energy = MinimizedEnergy(f'{build_path}/apo')
        lig_energy = MinimizedEnergy(f'{build_path}/lig')
        com_energy = MinimizedEnergy(f'{build_path}/com')
        diff_energy = com_energy - lig_energy - rec_energy
    except:
        success = False
    
    if not three_traj:
        print("1-traj calculation not ready")

    with open(f'{outpath}/metrics.csv','r') as metrics:
        dat = metrics.readlines()
    with open(f'{outpath}/metrics.csv','w') as metrics:
        metrics.write(dat[0].replace('\n',',Minimize,Minimize_U\n'))
        if success:
            metrics.write(dat[1].replace('\n',',{},{}\n'.format(diff_energy,0)))
        else:
            metrics.write(dat[1].replace('\n',',NA,NA\n'))


def RunMMGBSA(path, niter=1000):
    import mmgbsa
    crds = {'lig':f'{path}/lig.inpcrd','apo':f'{path}/apo.inpcrd','com':f'{path}/com.inpcrd'}
    prms = {'lig':f'{path}/lig.prmtop','apo':f'{path}/apo.prmtop','com':f'{path}/com.prmtop'}
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

def GetMetrics(smiles, method):
    """Method can be
        'dock'
        'minimize'
        'mmgbsa'
        'absolute' - Not in use
    """
    if method == 'dock':
        RunDocking(smiles,'input/com_axitinib.pdb','build')
        ParameterizeSystem('build')
    if method == 'minimize':
        RunDocking(smiles,'input/com_axitinib.pdb','build')
        ParameterizeSystem('build')
        RunMinimization('build','build')
    if method == 'mmgbsa':
        RunDocking(smiles,'input/com_axitinib.pdb','build')
        ParameterizeSystem('build')
        RunMinimization('build','build')
        RunMMGBSA('build', 500)
