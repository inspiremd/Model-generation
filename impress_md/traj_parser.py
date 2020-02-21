import sys, os
from contextlib import contextmanager

@contextmanager
def extractconfs(inpath, outpath, comp):
    """
    Extracts coordinates from a trajectory file (written by OpenMM in PDB format) and saves them as separate pdb files.
    """
    inpdb = inpath + '/' + comp + '_system_nosol.pdb'
    with open(inpdb,'r') as traj:
        lines = traj.readlines()
        cnt = 0
    
    for line in lines:
        if line.split()[0] in ('REMARK', 'END'):
            continue
        elif line.split()[0] == 'MODEL':
            cnt += 1
            outpdb = outpath + '/' + comp + '_MMGBSA_snap' + str(cnt) + '.pdb'
            f = open(outpdb, 'w')
            f.write(line)
        elif line.split()[0] ==  'ENDMDL':
            f.write(line)
            f.write("END")
            f.close()
        else:
            f.write(line)
    
def splitpdb(file1, file2, file3):
    """
    Takes a PDB file as input and splits it into two with one containing only the last residue while the other containing the rest.
    file1 is the input pdb. The last residue is written out to file3 and the rest is written out to file2.
    """    
    lnums = []

    with open(file1, 'r') as f1:
        for num, line in enumerate(f1, 1):
            if 'TER ' in line:
                lnums.append(num)

    cutoff = lnums[-2]
    f2 = open(file2, 'w')
    f3 = open(file3, 'w')

    with open(file1, 'r') as f1:
        for num, line in enumerate(f1):
            if num < cutoff:
                f2.write(line)
            else:
                f3.write(line)

    f2.write("END")
    f2.close()
    f3.close()


