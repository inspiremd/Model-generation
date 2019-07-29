import sys, os
sys.path.append('/home/adriandlee/IMPRESS/workflow/scripts')
import conf_gen, dock_conf
from openeye import oechem, oedocking
"""
Input
Output
Notes
"""

BuildFolders = True # this is just a workflow question

outpath = 'output/metrics.csv'
inpath  = 'input/smiles.txt'
molpath = 'build/ligands'
ifs = oechem.oemolistream()
ofs = oechem.oemolostream()


### Read in a text file that is <ID> <SMILEs>
#   ID could be a name or a number
def getSmiles(filename):
    smiles = {}
    with open(filename) as f:
        for line in f:
            (key,val) = line.split()
            smiles[key] = val
    return smiles

def selectEnantiomer(mol_list):
    return mol_list[0]

smiles = getSmiles(inpath)
mols = {}
for key in smiles:
    print(key,"has",len(conf_gen.fromString(smiles[key])),"enantiomers")
    mols[key] = selectEnantiomer(conf_gen.fromString(smiles[key]))
    mols[key].SetTitle(key)

receptor = oechem.OEGraphMol()
oedocking.OEReadReceptorFile(receptor,"build/receptor.oeb")

ligs = {}
docs = {}
for key in mols:
    print("On key",key)
    docs[key], ligs[key] = dock_conf.dockConf(receptor,mols[key])
    


# Write the best ligand to mol2 files
if not os.path.exists(molpath):
    os.mkdir(molpath)

for key in mols:
    if BuildFolders:
        lig_path = f'{molpath}/{key}/{key}.pdb'
        os.mkdir(f'{molpath}/{key}')
    else:
        lig_path = f'{molpath}/{key}.pdb'
    
    if ofs.open(lig_path):
        tmp_mol = list(ligs[key].GetConfs())[0]
        oechem.OEWriteMolecule(ofs, tmp_mol)

# Now analyze 
# Saves to a csv
import csv
# If no output file exists, make one
if not os.path.exists(outpath):
    with open(outpath,'w+') as ofs:
        ofs.write("ID,smiles,score_best,score_5\n")

for key in mols:
    line = dock_conf.getMetrics(docs[key],ligs[key])
    with open(outpath,'a') as ofs:
        wr = csv.writer(ofs)
        wr.writerow(line)
# Saves to a PANDAS df instead? Then adds to that after minimizing
