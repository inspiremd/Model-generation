import sys
from openeye import oechem, oeomega

def FromMol(mol, isomer=True):
    """
    Generates a set of conformers as an OEMol object
    Inputs:
        mol is an OEMol
        isomers is a boolean controling whether or not the various diasteriomers of a molecule are created
    """
    omegaOpts = oeomega.OEOmegaOptions()
    omega = oeomega.OEOmega(omegaOpts)
    out_conf = []
    
    if not isomer:
        ret_code = omega.Build(mol)
        if ret_code == oeomega.OEOmegaReturnCode_Success:
            out_conf.append(mol)
        else:
            oechem.OEThrow.Warning("%s: %s" % (mol.GetTitle(), oeomega.OEGetOmegaError(ret_code)))
    
    elif isomer:
        for enantiomer in oeomega.OEFlipper(mol.GetActive(),12,True):
            enantiomer = oechem.OEMol(enantiomer)
            ret_code = omega.Build(enantiomer)
            if ret_code == oeomega.OEOmegaReturnCode_Success:
                out_conf.append(enantiomer)
            else:
                oechem.OEThrow.Warning("%s: %s" % (mol.GetTitle(), oeomega.OEGetOmegaError(ret_code)))

    return out_conf

def FromString(smiles, isomer=True):
    """
    Generates an set of conformers from a SMILES string
    """
    mol = oechem.OEMol()
    if not oechem.OESmilesToMol(mol,smiles):
        print("SMILES invalid for string", smiles)
        return None
    else:
        return FromMol(mol,isomer)

# This is clearly not the right way to do it. TODO
def SelectEnantiomer(mol_list):
    return mol_list[0]
