# This script performs one pipeline task 

# SMILES should come from a database or from the caller of this pipeline.
SMILES="OCCCNc1cc(-c2ncnc(Nc3cccc(Cl)c3)n2)ccn1"

#query a python callable/server here to see if this step should be done or exit
time python docking.py -o "test/" -s $SMILES -i "input/"

#query a python callable/server here to see if this step should be done or exit
time python param.py -i "test"
time python mmgbsa.py -p "test" -n 0

#query a python callable/server here to see if this step should be done or exit
# This step should also ask for $n$ value.
time python mmgbsa.py -p "test" -n 1

#etc, we may add another step here. 
