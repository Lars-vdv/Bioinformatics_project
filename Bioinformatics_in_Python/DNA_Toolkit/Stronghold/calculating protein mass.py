import os
print(os.getcwd())  # Check current working directory

from Bioinformatics_in_Python.DNAToolkit import *



protein_seq = open(r"Bioinformatics_in_python\Stronghold\test_data\rosalind_prtm.txt",'r')
type(protein_seq)
protein_seq = protein_seq.read().replace('\n','')  # Read the file  
print(protein_seq)

def MW_protein_seq(protein_seq):
    """Calculates the molecular weight of a protein sequence"""
    return sum(MW_proteins[AA] for AA in protein_seq)


print(round(MW_protein_seq(protein_seq),ndigits=3))
