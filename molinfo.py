from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdchem
import numpy as np

class molinfo:
   
    # Constructor
    def __init__(self,mol):
        self.mol = mol
        self.allylic_carbons = []
        self.vinylic_carbons = []
        self.vinylic_pairs = set()

        self.primary_carbons_sp2 = []
        self.secondary_carbons_sp2 = []
        self.tertiary_carbons_sp2 = []
        self.most_to_least_sp2 = []

        self.primary_carbons_sp3 = []
        self.secondary_carbons_sp3 = []
        self.tertiary_carbons_sp3 = []
        self.most_to_least_sp3 = []

        self.collectinfo()
        self.find_most_to_least_sub_sp2()
        self.find_most_to_least_sub_sp3()
        
        
        
    def collectinfo(self):

        for atom in self.mol.GetAtoms():

            #check for tertiary/secondary/primary carbons
            if atom.GetAtomicNum() == 6:  #and str(atom.GetHybridization()) == "SP3" : there can be sp2 hy
                carbon_counter = 0
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        carbon_counter += 1

                if str(atom.GetHybridization()) == "SP3":
                    if carbon_counter == 1:
                        self.primary_carbons_sp3.append(atom.GetIdx())

                    elif carbon_counter == 2: 
                        self.secondary_carbons_sp3.append(atom.GetIdx())

                    elif carbon_counter == 3:
                            self.tertiary_carbons_sp3.append(atom.GetIdx())
                
                elif str(atom.GetHybridization()) == "SP2":
                    if carbon_counter == 1:
                        self.primary_carbons_sp2.append(atom.GetIdx())

                    elif carbon_counter == 2: 
                        self.secondary_carbons_sp2.append(atom.GetIdx())

                    elif carbon_counter == 3:
                            self.tertiary_carbons_sp2.append(atom.GetIdx())

        
            #check for vinylic carbons 
            if atom.GetAtomicNum() == 6 and str(atom.GetHybridization()) == "SP2":  #target atom is carbon and sp3
                for neighbor in atom.GetNeighbors():
                    bond = self.mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())

                    #check if neighbor is a carbon, if the bond type is double
                    if neighbor.GetAtomicNum() == 6 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        vinylic_pair = (min(atom.GetIdx(),neighbor.GetIdx()), max(atom.GetIdx(),neighbor.GetIdx()))
                        self.vinylic_pairs.add(vinylic_pair)
                        self.vinylic_carbons.append(atom.GetIdx())
                        #self.vinylic_carbons.append(neighbor.GetIdx())

                        for neighbor2 in atom.GetNeighbors():
                            
                            if neighbor2.GetAtomicNum() == 6 and str(neighbor2.GetHybridization()) == "SP3":

                                self.allylic_carbons.append(neighbor2.GetIdx())

    #Getters
    def get_allylic_carbons(self):

        return self.allylic_carbons

    def get_vinylic_carbons(self):

        return self.vinylic_carbons

    def get_vinylic_pairs(self):
        return self.vinylic_pairs

    def get_tertiary_carbons(self):

        return self.tertiary_carbons

    def get_secondary_carbons(self):

        return self.secondary_carbons

    def get_primary_carbons(self):

        return self.primary_carbons

    def get_most_to_least_sub_sp2(self):
        
        return self.most_to_least_sp2

    def get_most_to_least_sub_sp3(self):
        
        return self.most_to_least_sp3


    #Finders
    def find_most_to_least_sub_sp2(self):
        self.most_to_least_sp2 = self.tertiary_carbons_sp2 + self.secondary_carbons_sp2 + self.primary_carbons_sp2

    #Finders
    def find_most_to_least_sub_sp3(self):
        self.most_to_least_sp3 = self.tertiary_carbons_sp3 + self.secondary_carbons_sp3 + self.primary_carbons_sp3




    

                        

        
