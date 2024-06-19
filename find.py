from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdchem
import random
import molecule
from apply import *

def most_to_least_sub_in_vinylic(mol, vinylic_pair):
    #retuns a tuple: first index is most substituted, second index is least substituted
    atom_idx1 = vinylic_pair[0]
    atom_idx2 = vinylic_pair[1]
    
    atom1 = mol.GetAtomWithIdx(atom_idx1)
    atom2 = mol.GetAtomWithIdx(atom_idx2)

    atom1_hybrid = convert_hybrid_to_int(str(atom1.GetHybridization()))
    atom2_hybrid = convert_hybrid_to_int(str(atom2.GetHybridization()))

    if atom1_hybrid < atom2_hybrid:
        return [atom_idx2, atom_idx1]
    else: 
        return [atom_idx1, atom_idx2]




def convert_hybrid_to_int(hybrid):
    match hybrid:
        case "SP3": return 3
        case "SP2": return 2
        case "SP": return 1
        case _: pass


def allylic_carbons(mol):
    allylic_carbons = []
    

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and str(atom.GetHybridization()) == "SP3":  #target atom is carbon and sp3
            for neighbor in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())

                #check if neighbor is a carbon, if the bond type is double
                if neighbor.GetAtomicNum() == 6 and molecule.has_db(neighbor) :
                    
                   
                    allylic_carbons.append(atom.GetIdx())
                    break

    return allylic_carbons
    #allternative
    # # Define a SMARTS pattern for an allylic carbon
    # allylic_carbon_smarts = '[CH2,CH3]-[CH]=[CH]'

    # # Create a molecule from the SMARTS pattern
    # allylic_carbon_mol = Chem.MolFromSmarts(allylic_carbon_smarts)

    # # Find all matches of the SMARTS pattern in the molecule
    # return mol.GetSubstructMatches(allylic_carbon_mol)


