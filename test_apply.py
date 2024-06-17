

import unittest
from molecule import *
import find
import apply


class TestMol(unittest.TestCase):

    @classmethod
    def setUp(self):
        #------------Input
        #smiles = "CC(C)(C)OC(=O)C1=C(C=C(C=C1)N)OC2=C(C=C(C=C2)N)C3=NC(=NC(=N3)N)N"  #myster mol

        #smiles = "CCO"
        
        smiles = "CC=CCC=CCCCC=CC"
        self.molecule = Chem.MolFromSmiles(smiles)

        # Print canonical SMILES
        print("canonical SMILES " + Chem.MolToSmiles(self.molecule))
    
    @classmethod
    def tearDown(self):
        pass

    def test_allylic_bromination(self):
        new_mol = apply.allylic_bromination(self.molecule)

        img = Draw.MolToImage(new_mol)
        # Display the image (works well in Jupyter notebooks)
        img.show()

        


if __name__ == '__main__':
    unittest.main()





