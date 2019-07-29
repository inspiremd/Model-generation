from openeye import oechem, oedocking

### Read in a text file that is <Lig> <SMILEs>
def dockConf(receptor, mol, MAX_POSES = 5):
    dock = oedocking.OEDock()
    dock.Initialize(receptor)
    lig = oechem.OEMol()
    err = dock.DockMultiConformerMolecule(lig,mol,MAX_POSES)
    print("Error?",err)
    return dock, lig

### Returns an array of length MAX_POSES from above. This is the range of scores
def ligandScores(dock, lig):
    return [ dock.ScoreLigand(conf) for conf in lig.GetConfs() ]

def scoreRange(dock,lig):
    tmp = ligandScores(dock,lig)
    return tmp[0],tmp[-1]

def getMetrics(dock,lig):
    """ Outputs a single line for a matrix
    name, smiles, score_max, score_min 
    """
    out = []
    out.append(lig.GetTitle())
    out.append(oechem.OEMolToSmiles(lig))
    out += list(scoreRange(dock,lig))
    # Add energy-minimized MMPBSA
    return out
