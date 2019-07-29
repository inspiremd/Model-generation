#!/bin/bash

# Paths
curr=`pwd`
pLig=build/ligands

## Make a better-named mol2
## Make a $mol.prep
## Make a $mol.frcmod
## prepare the tleap.in


# First prepare receptor
## May need to fix charges
cd $pLig/..
echo "" > tleap.in
echo "source leaprc.protein.ff14SBonlysc" >> tleap.in
echo "source leaprc.gaff" >> tleap.in
echo "set default PBRadii mbondi3" >> tleap.in
echo "rec = loadPDB apo.pdb" >> tleap.in
echo "saveAmberParm rec apo.prmtop apo.inpcrd" >> tleap.in
echo "quit" >> tleap.in
tleap -f tleap.in

# For each file in $pLIG, run antechamber & tleap
for f in $pLig/*; do
    cd $f
    mol=${f#$pLig/}
    echo "Parameterizing $mol in $f"
    antechamber -i ${mol}.pdb -fi pdb -o ${mol}.mol2 -fo mol2 -c bcc -pf y
    parmchk -i ${mol}.mol2 -f mol2 -o ${mol}.frcmod
    source $main/scripts/prepare_leap.sh $mol
    tleap -f tleap.in
done

cd $curr
