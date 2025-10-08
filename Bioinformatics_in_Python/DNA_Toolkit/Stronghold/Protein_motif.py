import os
print(os.getcwd())  # Check current working directory

from bio_seq import *
from bio_structs import *
from utilities import *
import urllib.request

seq_1 = urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/B5ZC00.fasta").read().decode('utf-8') # Fetch sequence data from URL, read and decode it to string
print(type(seq_1))
print(seq_1)

with open("test_data/rosalind_mprt.txt", "r") as f:
    temp = f.readlines()
print(temp)



for i in file: 
    
    with open("test_data/rosalind_mprt.txt", "r") as f:
        temp = f.readlines()

