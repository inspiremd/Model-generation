from openeye import oechem, oedocking

def PrepareReceptor(pdb,padding=4,outpath=""):
    """
    Prepares a receptor from a pdb with a crystalized ligand
    Padding controls the docking region.
    If outpath is given, PrepareReceptor will write an openeye binary (oeb) of the receptor structure. This will be faster than rebuilding the receptor every time.
    """
    com = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    if ifs.open(pdb):
        oechem.OEReadPDBFile(ifs, com)
        ifs.close()
    oechem.OEDetermineConnectivity(com)
    nparts, connect = oechem.OEDetermineComponents(com)
    if(nparts != 2):
        print("ERR in dock_conf::prepareReceptor. PDB doesn't have 2 connected components")
        ## TODO: What is a good way to catch errors?
    # Get apo
    pred = oechem.OEPartPredAtom(connect)
    pred.SelectPart(nparts)
    lig = oechem.OEGraphMol()
    oechem.OESubsetMol(lig, com, pred)
    
    # Get protein
    pred = oechem.OEPartPredAtom(connect)
    pred.SelectPart(1)
    prot = oechem.OEGraphMol()
    oechem.OESubsetMol(prot, com, pred)
    
    # Get box dimensions by iterating over ligand
    x_min = y_min = z_min = float('inf')
    x_max = y_max = z_max = -float('inf')
    crd = lig.GetCoords()
    for atm in crd:
        x,y,z = crd[atm]
        if x < x_min:
            x_min = x
        if y < y_min:
            y_min = y
        if z < z_min:
            z_min = z
        if x > x_max:
            x_max = x
        if y > y_max:
            y_max = y
        if z > z_max:
            z_max = z
    x_min -= padding
    y_min -= padding
    z_min -= padding
    x_max += padding
    y_max += padding
    z_max += padding
    
    # Now prepare the receptor
    receptor = oechem.OEGraphMol()
    box = oedocking.OEBox(x_max, y_max, z_max, x_min, y_min, z_min)
    oedocking.OEMakeReceptor(receptor, prot, box)
    
    if not outpath == "":
        oedocking.OEWriteReceptorFile(receptor,f'{outpath}/receptor.oeb')
    return receptor

def PrepareReceptorFromBinary(filename):
    receptor = oechem.OEGraphMol()
    oedocking.OEReadReceptorFile(receptor,filename)
    return receptor

def DockConf(receptor, mol, MAX_POSES = 5):
    dock = oedocking.OEDock()
    dock.Initialize(receptor)
    lig = oechem.OEMol()
    err = dock.DockMultiConformerMolecule(lig,mol,MAX_POSES)
    # print("Error?",err)
    return dock, lig

def WriteStructures(receptor, lig, apo_path, lig_path):
    ofs = oechem.oemolostream()
    success = True
    if ofs.open(apo_path):
        oechem.OEWriteMolecule(ofs,receptor)
        ofs.close()
    else:
        success = False
    # If MAX_POSES != 1, we should select the top pose to save
    conf = list(lig.GetConfs())[0]
    if ofs.open(lig_path):
        oechem.OEWriteMolecule(ofs,conf)
        ofs.close()
    else:
        success = False
    return success

### Returns an array of length MAX_POSES from above. This is the range of scores
def LigandScores(dock, lig):
    return [ dock.ScoreLigand(conf) for conf in lig.GetConfs() ]

def BestDockScore(dock, lig):
    return LigandScores(dock,lig)[0]

def ScoreRange(dock,lig):
    tmp = LigandScores(dock,lig)
    return tmp[0],tmp[-1]

### Not in use
def GetMetrics(dock,lig):
    """ Outputs a single line for a matrix
    name, smiles, score_max, score_min 
    """
    out = []
    out.append(lig.GetTitle())
    out.append(oechem.OEMolToSmiles(lig))
    out += list(ScoreRange(dock,lig))
    # Add energy-minimized MMPBSA
    return out
