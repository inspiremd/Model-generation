from openeye import oechem, oedepict, oeomega

"""
Example script that generates conformers for a small set of molecules
"""

smiles = {'cbd':'CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O',
            'cbc':'CCCCCC1=CC2=C(C=CC(O2)(C)CCC=C(C)C)C(=C1)O',
            'cbn':'CCCCCC1=CC2=C(C(=C1)O)C3=C(C=CC(=C3)C)C(O2)(C)C',
            'thc':'CCCCCC1=CC2=C(C3C=C(CCC3C(O2)(C)C)C)C(=C1)O',
            'cbg':'CCCCCC1=CC(=C(C(=C1)O)CC=C(C)CCC=C(C)C)O'}

### Make the basic molecules from the smiles strings
mols = dict()
for mol in smiles:
    mols[mol] = oechem.OEMol()
    if not oechem.OESmilesToMol(mols[mol],smiles[mol]) : print("SMILES invalid for mol",mol)


import conf_gen
confs = dict()
for mol in mols:
    confs[mol] = conf_gen.fromMol(mols[mol])
    # confs is now a dictionary of lists. 
    # The length of each list is the number of isomers

### Print first 5 conformers of each molecule
ofs = oechem.oemolostream()
for mol in mols:
    for isomer in enumerate(confs[mol]): # Isomer number
        for conf in enumerate(isomer[1].GetConfs()):
            ## Normally these are treated as a single object, not separate conformers
            if (conf[0] < 5):
                if ofs.open(mol+"-isomer_"+str(isomer[0])+"-conf_"+str(conf[0])+".mol2"):
                    oechem.OEWriteMolecule(ofs,conf[1])
                    ofs.close()
