from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdchem
import random
import molecule
from apply import *

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

