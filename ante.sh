#!/bin/bash
start=`date +%s`

path=$1
prev=`pwd`
cd $path
antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y
parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod

# Wrap tleap to get $path/(com|lig|apo).(inpcrd|prmtop)
echo "source leaprc.protein.ff14SBonlysc" > leap.in
echo "source leaprc.gaff" >> leap.in
echo "set default PBRadii mbondi3" >> leap.in
echo "rec = loadPDB apo.pdb # May need full filepath?" >> leap.in
echo "saveAmberParm rec apo.prmtop apo.inpcrd" >> leap.in
echo "lig = loadmol2 lig.mol2" >> leap.in
echo "loadAmberParams lig.frcmod" >> leap.in
echo "com = combine {rec lig}" >> leap.in
echo "saveAmberParm lig lig.prmtop lig.inpcrd" >> leap.in
echo "saveAmberParm com com.prmtop com.inpcrd" >> leap.in
echo "quit" >> leap.in

tleap -f leap.in
end=`date +%s`
runtime=$((end-start))

echo "Execution time: $runtime" > ante.log

cd $prev
