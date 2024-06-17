from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdchem
import random
import molecule
import find

def allylic_bromination(mol):
    # Apply bromination to the specified allylic carbon
    allylic_carbons = find.allylic_carbons(mol)

    # Randomly pick one allylic carbon index
    atom_to_modify = random.choice(allylic_carbons)

    molecule.getAtomInfo(mol, atom_to_modify)

    return brominate_atom(mol,atom_to_modify)

def brominate_atom(molecule, atom_idx):
    """
    Brominate the given atom in the molecule by adding a bromine (Br) atom.
    
    Parameters:
    molecule (rdkit.Chem.rdchem.RWMol): The molecule to be modified.
    atom_idx (int): The index of the atom to brominate.
    
    Returns:
    rdkit.Chem.rdchem.RWMol: The modified molecule.
    """
    # Convert to RWMol (a modifiable molecule)
    rw_mol = Chem.RWMol(molecule)
    
    # Create a new bromine atom
    br_atom = Chem.Atom(35)  # Atomic number for Bromine is 35
    
    # Add the bromine atom to the molecule
    br_idx = rw_mol.AddAtom(br_atom)
    
    # Create a bond between the existing atom and the new bromine atom
    rw_mol.AddBond(atom_idx, br_idx, rdchem.BondType.SINGLE)
    
    # Return the modified molecule
    return rw_mol
