from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdchem
import random
import molecule
import find
from molinfo import molinfo

import collections

def allylic_bromination(mol, info):
    # Apply bromination to the specified allylic carbon
    allylic_carbons = info.get_allylic_carbons()

    # Randomly pick one allylic carbon index
    atom_to_modify = random.choice(allylic_carbons)

    molecule.getAtomInfo(mol, atom_to_modify)

    # Convert to RWMol (a modifiable molecule)
    rw_mol = Chem.RWMol(mol)

    return brominate_atom(rw_mol,atom_to_modify)

def hydrochlorination(mol,info):
    vinylic_carbons = info.get_vinylic_carbons()

    most_to_least = info.get_most_to_least_sub()

    intersection = list(collections.Counter(vinylic_carbons) & collections.Counter(most_to_least))

    print(vinylic_carbons)
    print(most_to_least)
    atom_idx1 = intersection[0]
    atom = mol.GetAtomWithIdx(atom_idx1)
    atom_idx2 = 0

    for bond in atom.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            # Get the other atom index connected by the double bond
            atom_idx2 = bond.GetOtherAtomIdx(atom_idx1)

    # Convert to RWMol (a modifiable molecule)
    rw_mol = Chem.RWMol(mol)

    rw_mol = chlorinate_atom(rw_mol, atom_idx1)

    rw_mol = break_pi_bond(rw_mol,atom_idx1, atom_idx2)

    return rw_mol

def break_pi_bond(rw_mol,atom_idx1, atom_idx2):

    # Change the bond order from double to single
    rw_mol.RemoveBond(atom_idx1, atom_idx2)
    rw_mol.AddBond(atom_idx1, atom_idx2, Chem.rdchem.BondType.SINGLE)

    # Adjust the hydrogen counts on the atoms to maintain valence
    atom1 = rw_mol.GetAtomWithIdx(atom_idx1)
    atom2 = rw_mol.GetAtomWithIdx(atom_idx2)
    #atom1.SetNumExplicitHs(atom1.GetTotalNumHs() + 1)
    #atom2.SetNumExplicitHs(atom2.GetTotalNumHs() + 1)

    # Update the molecule structure
    Chem.SanitizeMol(rw_mol)
    return rw_mol

def chlorinate_atom(rw_mol,atom_idx):
 
   
    # Create a new bromine atom
    cl_atom = Chem.Atom(17)  # Atomic number for Chlorine is 17
    
    # Add the bromine atom to the molecule
    cl_idx = rw_mol.AddAtom(cl_atom)
    
    # Create a bond between the existing atom and the new bromine atom
    rw_mol.AddBond(atom_idx, cl_idx, rdchem.BondType.SINGLE)
    
    # Return the modified molecule
    return rw_mol

def brominate_atom(rw_mol, atom_idx):
    """
    Brominate the given atom in the molecule by adding a bromine (Br) atom.
    
    Parameters:
    molecule (rdkit.Chem.rdchem.RWMol): The molecule to be modified.
    atom_idx (int): The index of the atom to brominate.
    
    Returns:
    rdkit.Chem.rdchem.RWMol: The modified molecule.
    """
    
    
    # Create a new bromine atom
    br_atom = Chem.Atom(35)  # Atomic number for Bromine is 35
    
    # Add the bromine atom to the molecule
    br_idx = rw_mol.AddAtom(br_atom)
    
    # Create a bond between the existing atom and the new bromine atom
    rw_mol.AddBond(atom_idx, br_idx, rdchem.BondType.SINGLE)
    
    # Return the modified molecule
    return rw_mol
