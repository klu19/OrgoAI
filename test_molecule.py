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
        
        #smiles = "CC=CCC=CCCCC=CC"
        smiles = "CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4=C3C=CC=C4)C"
        self.molecule = Chem.MolFromSmiles(smiles)

        # Print canonical SMILES
        print("canonical SMILES " + Chem.MolToSmiles(self.molecule))
    
    @classmethod
    def tearDown(self):
        pass

    #@unittest.skip("Skipping find allylic carbons test")
    def test_Draw(self):
        # Draw the molecule
        img = Draw.MolToImage(self.molecule)

        # Display the image (works well in Jupyter notebooks)
        img.show()
    def test_tert(self):
        
        print(molecule.getTertiary(self.molecule))
        

    


if __name__ == '__main__':
    unittest.main()





