import sys, os
sys.path.append('/home/adriandlee/IMPRESS/workflow/scripts')
import minimize

rec = 'build/apo'
inpath = 'build/ligands'
outpath = 'output/metrics.csv'

mols = []
with open('input/smiles.txt') as ifs:
    for line in ifs:
        (key,val) = line.split()
        mols.append(key)

### Get energies
rec_energy = minimize.MinimizedEnergy(rec)
lig_energy = {}
com_energy = {}
diff_energy = {}
for key in mols:
    try:
        lig_energy[key] = minimize.MinimizedEnergy(f'{inpath}/{key}/ligand')
        com_energy[key] = minimize.MinimizedEnergy(f'{inpath}/{key}/complex')
        diff_energy[key] = com_energy[key] - lig_energy[key] - rec_energy
    except:
        continue

# Append difference to column in metrics.csv
import csv
with open(outpath,'r') as csv_in:
    with open(f'{outpath}.temp','w') as csv_out:
        writer = csv.writer(csv_out,lineterminator='\n')
        reader = csv.reader(csv_in)
        dat = []
        
        row = next(reader)
        row.append('Minimized Energy')
        dat.append(row)
        
        for row in reader:
            try:
                row.append(diff_energy[row[0]])
            except:
                row.append('NA')
            dat.append(row)
        writer.writerows(dat)
        # assumes that row[0] is the key
        # Returns 'NA' when the key is not found

