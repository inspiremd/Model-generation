import argparse

def RunDocking_(smiles, inpath, outpath, padding=4):
    from . import conf_gen
    from . import dock_conf
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    confs = conf_gen.SelectEnantiomer(conf_gen.FromString(smiles))
    # This receptor can be pre-compiled to an oeb. It speeds things up
    filename, file_extension = os.path.splitext(inpath)
    #if file_extension == ".oeb":
    #    receptor = dock_conf.PrepareReceptorFromBinary(inpath)
    #else: # else it is a pdb
    #    receptor = dock_conf.PrepareReceptor(inpath,padding,outpath)

    dock, lig, receptor = dock_conf.DockConf("input/receptor.oeb",confs,MAX_POSES=1)

    # Currently we generate 200 conformers for each ligand, but only take
    #   the best pose, as scored by Openeye. It may be useful to consider
    #   something about the range of poses.

    with open(f'{outpath}/metrics.csv','w+') as metrics:
        metrics.write("Dock,Dock_U\n")
        metrics.write("{},{}\n".format(dock_conf.BestDockScore(dock,lig),0))
    dock_conf.WriteStructures(receptor, lig, f'{outpath}/apo.pdb', f'{outpath}/lig.pdb')
    # # If you uncomment the three lines below, it will save an image of the 2D
    #   molecule. This is useful as a sanity check.
    # from openeye import oedepict
    # oedepict.OEPrepareDepiction(lig)
    # oedepict.OERenderMolecule(f'{outpath}/lig.png',lig)
    return dock_conf.BestDockScore(dock,lig)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smile", type=str, default=None, required=True)
    parser.add_argument("--inpath", type=str, defualt=None, required=True)
    parser.add_argument('--outpath', type=str, default=None, required=True)

if __name__ == '__main__':
    args = get_args()
    RunDocking_(args.smiles, args.inpath, args.outpath)