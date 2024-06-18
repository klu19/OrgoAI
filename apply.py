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

    return add_atom(rw_mol,atom_to_modify,35,1)

def radical_halogenation(mol, info, halo):

    most_to_least = info.get_descending_rad()
    print(most_to_least)
    #assumes there is a reactive site
    atom_idx = most_to_least[0]
    
    # Convert to RWMol (a modifiable molecule)
    rw_mol = Chem.RWMol(mol)

    rw_mol = add_atom(rw_mol, atom_idx, halo, 1) #halogenate

    return rw_mol

def hydrohalogenation(mol,info,halo):
    vinylic_carbons = info.get_vinylic_carbons()

    most_to_least = info.get_most_to_least_sub_sp2()

    intersection = list(collections.Counter(vinylic_carbons) & collections.Counter(most_to_least))

    print(vinylic_carbons)
    print(info.get_vinylic_pairs())
    print(most_to_least)
    atom_idx1 = intersection[0]
    atom = mol.GetAtomWithIdx(atom_idx1)
    atom_idx2 = 0

    #scuffed, just use vinylic pair
    for bond in atom.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            # Get the other atom index connected by the double bond
            atom_idx2 = bond.GetOtherAtomIdx(atom_idx1)

    # Convert to RWMol (a modifiable molecule)
    rw_mol = Chem.RWMol(mol)

    rw_mol = add_atom(rw_mol, atom_idx1,halo, 1)

    rw_mol = break_pi_bond(rw_mol,atom_idx1, atom_idx2)

    return rw_mol

def dihalogenation(mol,info,halo):

    vinylic_pairs = info.get_vinylic_pairs()

    # Convert set to list (since sets do not support indexing)
    pairs_list = list(vinylic_pairs)

    # Randomly pick one vinylic pair index
    pair_to_modify = random.choice(pairs_list)

    atom_idx1 = pair_to_modify[0]

    atom_idx2 = pair_to_modify[1]
    # Convert to RWMol (a modifiable molecule)
    rw_mol = Chem.RWMol(mol)
    rw_mol = add_atom(rw_mol, atom_idx1, halo, 1) #chlorinate
    rw_mol = add_atom(rw_mol, atom_idx2, halo,1) #chlorinate

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

#bond direction, else is obsolete


#     bond = rw_mol.GetBondBetweenAtoms(atom_idx, cl_idx)

#     match bondir: 
#         case 1 :
#             # Modify the bond to be drawn as a wedge (up) or dash (down)
#             print("enntered this one")
#             bond.SetBondDir(Chem.rdchem.BondDir.ENDUPRIGHT) # Wedge (up) direction
#         case 2 :

#             print("then this one")
#             bond.SetBondDir(Chem.rdchem.BondDir.ENDDOWNRIGHT)  # Dash (down) direction
#         case _:
#             pass


def add_atom(rw_mol, atom_idx, new_atomic_num, bond):

    new_atom = Chem.Atom(new_atomic_num)

    new_atom_idx = rw_mol.AddAtom(new_atom)

    # Create a bond between the existing atom and the new atom
    btype = rdchem.BondType.SINGLE
    match bond:
        case 1: btype = rdchem.BondType.SINGLE
        case 2: btype = rdchem.BondType.DOUBLE
        case 3: btype = rdchem.BondType.TRIPLE
        case _: pass

    rw_mol.AddBond(atom_idx, new_atom_idx, btype)

    return rw_mol