

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

        #smiles = "CC(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C=C)C(C)(C)C"

        #smiles = "C=CCC(CCC)CC1=CC=CC=C1"

        #smiles = "CC1=CC=CC=C1CC2=CC=CC=C2CC3=CC=CC=C3"

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
      
        
        new_mol = apply.hydrohalogenation(self.molecule, self.mol1,17)

        img = Draw.MolToImage(new_mol)
        # Display the image (works well in Jupyter notebooks)
        img.show()
    @unittest.skip("Skipping apply dichlorination test")
    def test_dichlorination(self):
      
        Draw.MolToImage(self.molecule).show()
        new_mol = apply.dihalogenation(self.molecule, self.mol1,17)
    
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
    
    @unittest.skip("Skipping radical halogenation test")
    def test_radical_bromination(self):
        
        new_mol = apply.radical_halogenation(self.molecule, self.mol1,35)
        Draw.MolToImage(self.molecule).show()
        Draw.MolToImage(new_mol).show()

    @unittest.skip("Skipping add hydroxyl test")
    def test_add_OH(self):
        rw_mol = Chem.RWMol(self.molecule)
        new_mol = apply.add_OH(rw_mol,0)
        Draw.MolToImage(self.molecule).show()
        Draw.MolToImage(new_mol).show()

    #@unittest.skip("Skipping haloh test")
    def test_halohydrin(self):
        
        rw_mol = Chem.RWMol(self.molecule)
        Draw.MolToImage(self.molecule).show()
        print(self.mol1.get_vinylic_carbons)
        new_mol = apply.halohydrin(rw_mol,self.mol1,17)
        
        Draw.MolToImage(new_mol).show()

if __name__ == '__main__':
    unittest.main()





