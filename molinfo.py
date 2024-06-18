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
        
        self.primary_carbons = []
        self.secondary_carbons = []
        self.tertiary_carbons = []
        self.most_to_least = []

        self.collectinfo()
        self.find_most_to_least_sub()
        
        
    def collectinfo(self):

        for atom in self.mol.GetAtoms():

            #check for tertiary/secondary/primary carbons
            if atom.GetAtomicNum() == 6:  #and str(atom.GetHybridization()) == "SP3" : there can be sp2 hy
                carbon_counter = 0
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        carbon_counter += 1

                if carbon_counter == 1:
                    self.primary_carbons.append(atom.GetIdx())

                elif carbon_counter == 2: 
                    self.secondary_carbons.append(atom.GetIdx())

                elif carbon_counter == 3:
                        self.tertiary_carbons.append(atom.GetIdx())

        
            #check for vinylic carbons 
            if atom.GetAtomicNum() == 6 and str(atom.GetHybridization()) == "SP2":  #target atom is carbon and sp3
                for neighbor in atom.GetNeighbors():
                    bond = self.mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())

                    #check if neighbor is a carbon, if the bond type is double
                    if neighbor.GetAtomicNum() == 6 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                
                        self.vinylic_carbons.append(atom.GetIdx())

                        for neighbor in atom.GetNeighbors():
                            
                            if neighbor.GetAtomicNum() == 6 and str(neighbor.GetHybridization()) == "SP3":

                                self.allylic_carbons.append(neighbor.GetIdx())

    #Getters
    def get_allylic_carbons(self):

        return self.allylic_carbons

    def get_vinylic_carbons(self):

        return self.vinylic_carbons

    def get_tertiary_carbons(self):

        return self.tertiary_carbons

    def get_secondary_carbons(self):

        return self.secondary_carbons

    def get_primary_carbons(self):

        return self.primary_carbons

    def get_most_to_least_sub(self):

        return self.most_to_least


    #Finders
    def find_most_to_least_sub(self):
        self.most_to_least = self.tertiary_carbons + self.secondary_carbons + self.primary_carbons



    

                        

        
