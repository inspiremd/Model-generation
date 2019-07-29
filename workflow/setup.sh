#!/bin/bash

if [ ! -d input ]; then
    echo "No input" 1>&2
fi

# Optional input is input complex

: '
Input for prepare_receptor is just a pdb named input/com.pdb OR user input. It only works when there are two chains, the former is the protein, and the latter is the ligand. It generates a few files in build (receptor.oeb, apo.pdb, and original_lig.pdb).

Input for prepare_dock is the receptor in build and input/smiles.txt. It outputs a build/ligands directory containing docked ligands at the right coordinates to pair with build/apo.*. 

Input for prepare_parm is the set of directory build/ligands/* and build/apo.pdb. It outputs (ligand|complex).(prmtop|inpcrd) into the dir of each respective ligand. apo.(prmtop|inpcrd) is also created in build/. These should be simulation ready?

Input for minimize_script are the prmtop + inpcrd files for ligand, complex, and apo receptor.
'

[[ -d output ]] || mkdir output

echo "Preparing receptor..."
[[ -e build/receptor.oeb ]] || python scripts/prepare_receptor.py $1
echo "Receptor built."

echo "Docking ligands..."
[[ -e output/metrics.csv ]] || python scripts/prepare_dock.py
echo "Ligands docked."

echo "Parameterizing ligands..."
./scripts/prepare_param.sh
echo "Receptor & ligands parameterized."

echo "Minimizing ligands..."
python scripts/minimize_script.py
mv output/metrics.csv.temp output/metrics.csv
echo "Minimized energy written"
