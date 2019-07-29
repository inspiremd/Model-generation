from openeye import oechem, oedocking
import sys, os
"""
Input
    Required input is p_complex, the pdb of a protein and ligand
    If no argument is given, assumes the path input/com.pdb
Output
    build/apo.pdb, build/original_lig.mol2, build/receptor.oeb
        apo.mol2 is the apo structure of the protein
        receptor.oeb is for the docking
Notes
    This won't work if the system has multiple chains
    The output apo protein doesn't display with NewCartoon in vmd.
    The padding around the existing ligand is arbitrarily set to 4 A. 2 A was too small
"""

# Inputs
if len(sys.argv) == 2:
    p_complex = sys.argv[1]
else:
    p_complex = "input/com.pdb"
pad = 4 # padding around ligand. I think 2 is too small


if not os.path.exists("build"):
    os.mkdir("build")

ifs = oechem.oemolistream()
ofs = oechem.oemolostream()

com = oechem.OEGraphMol()
if ifs.open(p_complex):
    oechem.OEReadPDBFile(ifs,com)
    ifs.close()

# Find the connected components of the system, remove the ligand
oechem.OEDetermineConnectivity(com)
nparts, connect = oechem.OEDetermineComponents(com)
print(f'Number of connected components: {nparts}')
## Raise error or warning if there are more than 2
##   which means you need to select the right ligand

# Separate the ligand
pred = oechem.OEPartPredAtom(connect)
pred.SelectPart(nparts)
lig = oechem.OEGraphMol()
oechem.OESubsetMol(lig, com, pred)
if ofs.open("build/original_lig.mol2"):
    oechem.OEWriteMolecule(ofs,lig)
    ofs.close()

# Assume there's only one protein chain,
del pred
pred = oechem.OEPartPredAtom(connect)
pred.SelectPart(1)
prot = oechem.OEGraphMol()
oechem.OESubsetMol(prot, com, pred)
if ofs.open("build/apo.pdb"):
    oechem.OEWriteMolecule(ofs,prot)
    ofs.close()


# Now iterate over these indices to create a box for docking
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

# Padding of 4 angstrom is not necessarily best. 2 A is too small
x_min -= pad
y_min -= pad
z_min -= pad
x_max += pad
y_max += pad
z_max += pad
print("Min: (%.2f, %.2f, %.2f)" % (x_min,y_min,z_min))
print("Max: (%.2f, %.2f, %.2f)" % (x_max,y_max,z_max))

# Now prepare the receptor
receptor = oechem.OEGraphMol()
box = oedocking.OEBox(x_max, y_max, z_max, x_min, y_min, z_min)
oedocking.OEMakeReceptor(receptor, prot, box)

# # Clean the receptor as needed
# if oedocking.OEReceptorHasBoundLigand(receptor):
#     lig = oedocking.OEReceptorGetBoundLigand(receptor)
#     ofs = oechem.oemolostream("original_lig.mol2")
#     oechem.OEWriteMolecule(ofs,lig)
#     ofs.close()
#     oedocking.OEReceptorClearBoundLigand(receptor)

oedocking.OEWriteReceptorFile(receptor,"build/receptor.oeb")

### So now we have a built receptor. 
# ISSUES:
#   The receptor doesn't recognize the ligand, neither as a ligand nor as an extra molecule
#   My choice of padding by 2 Angstrom is arbitrary.
#   This will break if the last connected component isn't the ligand
#   I don't know what would happen if the protein has multiple chains.








