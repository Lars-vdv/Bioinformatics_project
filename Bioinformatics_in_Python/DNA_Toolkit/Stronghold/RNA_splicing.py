import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

print(os.getcwd())
file_path = r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\test_data\rosalind_splc.txt"

from bio_seq import bio_seq
from bio_structs import *
from utilities import *
import re

gene_and_introns = bio_seq(seq="ACTG", seq_type="DNA")
gene_and_introns = parse_fasta(path = file_path, seq_type="DNA", create_variables=True)


v = list(gene_and_introns.values())

print(gene_and_introns["seq_1"])

#get the sequence of the objects in the dictionary
for i in range(len(gene_and_introns)):
    name = str(f"seq_{i+1}")
    
    gene_and_introns[name] = getattr(gene_and_introns[name],"seq")



 