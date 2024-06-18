
import unittest
from molecule import *
import find
import apply


class TestMolinfo(unittest.TestCase):

    @classmethod
    def setUp(self):
        #------------Input
        #smiles = "CC(C)(C)OC(=O)C1=C(C=C(C=C1)N)OC2=C(C=C(C=C2)N)C3=NC(=NC(=N3)N)N"  #myster mol

        #smiles = "CCO"
        
        #smiles = "CC=CCC=CCCCC=CC"

        #smiles = 'c1cc(C(=O)O)c(OC(=O)C)cc1'

        # Example molecule: Toluene
        #smiles = 'CC1=CC=CC=C1'

        smiles = 'CC(Cc1ccccc1)OC2=CC(COC(Cc3ccccc3)C)C(COC(Cc4ccccc4)C)=C2'

        self.molecule = Chem.MolFromSmiles(smiles)

        # Print canonical SMILES
        print("canonical SMILES " + Chem.MolToSmiles(self.molecule))
    
    @classmethod
    def tearDown(self):
        pass

    @unittest.skip("Skipping find allylic carbons test")
    def test_find_allylic_carbons(self):
        # allylic_carbons = find.allylic_carbons(self.molecule)
        # print("Allylic carbon indices:", ", ".join(map(str, allylic_carbons)))
        print(find.allylic_carbons(self.molecule))

    #@unittest.skip("Skipping draw test")
    def test_Draw(self):
        # Draw the molecule
        img = Draw.MolToImage(self.molecule)

        # Display the image (works well in Jupyter notebooks)
        img.show()

    def test_molinfo(self):
        mol1 = molinfo(self.molecule)
        

if __name__ == '__main__':
    unittest.main()





