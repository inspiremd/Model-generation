#!/bin/bash
### This script assumes that we're in the folder for the ligand build.
### $1 is the name of the ligand
#   Input is 
#   Output is a complex.inpcrd, complex.prmtop, ligand.inpcrd, ligand.prmtop

mol=$1

# Echo everything into a leap.in file in the current directory
echo "" > tleap.in
echo "source leaprc.protein.ff14SBonlysc" >> tleap.in
echo "source leaprc.gaff" >> tleap.in
echo "set default PBRadii mbondi3" >> tleap.in
echo "lig = loadmol2 ${mol}.mol2" >> tleap.in
echo "loadAmberParams ${mol}.frcmod" >> tleap.in
echo "rec = loadPDB ../../apo.pdb" >> tleap.in
echo "com = combine {rec lig}" >> tleap.in
echo "saveAmberParm lig ligand.prmtop ligand.inpcrd" >> tleap.in
echo "saveAmberParm com complex.prmtop complex.inpcrd" >> tleap.in
echo "quit" >> tleap.in