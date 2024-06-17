from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdchem
import random
from find import *
from apply import *


def molecule_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    # Convert RDKit molecule to graph representation
    # Each atom is a node, each bond is an edge
    return mol

def has_db(atom):
    """
    Check if the given atom has any double bonds.
    
    Parameters:
    atom (rdkit.Chem.rdchem.Atom): The atom to check.
    
    Returns:
    bool: True if the atom has at least one double bond, False otherwise.
    """
    for bond in atom.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            return True
    return False



def print_atom_bonds(molecule, atom_idx):

    atom = molecule.GetAtomWithIdx(atom_idx)
    print(f"Atom index {atom_idx} (Atomic number {atom.GetAtomicNum()}) bonds:")
    for bond in atom.GetBonds():
        bond_type = bond.GetBondType()
        other_atom_idx = bond.GetOtherAtomIdx(atom_idx)
        other_atom = molecule.GetAtomWithIdx(other_atom_idx)
        other_atom_num = other_atom.GetAtomicNum()
        print(f" - Bond type: {bond_type}, Connected to atom index: {other_atom_idx} (Atomic number {other_atom_num})")


def getAtomInfo(mol,atom_idx):
    atom = mol.GetAtomWithIdx(atom_idx)
    print("Atomic Number: " + str(atom.GetAtomicNum()))
    print("Hybridization: " + str(atom.GetHybridization()))
    print_atom_bonds(mol, atom_idx)
    

    




    
import heapq

def find_synthesis_path(start_smiles, target_smiles):
    start_mol = molecule_to_graph(start_smiles)
    target_mol = molecule_to_graph(target_smiles)
    
    # Priority queue for A* search
    pq = [(0, start_mol)]
    visited = set()
    
    while pq:
        cost, current_mol = heapq.heappop(pq)
        if current_mol in visited:
            continue
        visited.add(current_mol)
        
        if current_mol == target_mol:
            return True
        
        # Generate next states using reaction rules
        for next_mol in generate_next_states(current_mol):
            heapq.heappush(pq, (cost + 1, next_mol))
    
    return False