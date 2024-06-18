

import unittest
from molecule import *
import find
import apply
from molinfo import molinfo

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


class TestMol(unittest.TestCase):

    @classmethod
    def setUp(self):
        #------------Input
        #smiles = "CC(C)(C)OC(=O)C1=C(C=C(C=C1)N)OC2=C(C=C(C=C2)N)C3=NC(=NC(=N3)N)N"  #myster mol

        #smiles = "CCO"
        
        #smiles = "CC=CCC=CCCCC=CC"

        smiles = "CC(=C)C=CC=C"

        self.molecule = Chem.MolFromSmiles(smiles)

        self.mol1 = molinfo(self.molecule)

        # Print canonical SMILES
        print("canonical SMILES " + Chem.MolToSmiles(self.molecule))
    
    @classmethod
    def tearDown(self):
        pass

    @unittest.skip("Skipping apply allylic carbons test")
    def test_allylic_bromination(self):
        
        
        new_mol = apply.allylic_bromination(self.molecule, self.mol1)

        img = Draw.MolToImage(new_mol)
        # Display the image (works well in Jupyter notebooks)
        img.show()
    @unittest.skip("Skipping apply hydrochlorination test")
    def test_hydrochlorination(self):
      
        
        new_mol = apply.hydrochlorination(self.molecule, self.mol1)

        img = Draw.MolToImage(new_mol)
        # Display the image (works well in Jupyter notebooks)
        img.show()

    def test_dichlorination(self):
      
        
        new_mol = apply.dichlorination(self.molecule, self.mol1)
        Draw.MolToImage(self.molecule).show()
        Draw.MolToImage(new_mol).show()

        #not working
        # # Visualize the modified molecule with the specified bond direction using RDKit drawing tools
        # drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)
        # drawer.DrawMolecule(new_mol)
        # drawer.FinishDrawing()

        # # Get the SVG representation
        # svg = drawer.GetDrawingText()

        # # Output or display the SVG image
        # with open('molecule.svg', 'w') as f:
        #     f.write(svg)
        
       

        


if __name__ == '__main__':
    unittest.main()





