import os
print(os.getcwd())

import re
import inspect



from Bioinformatics_in_Python.DNAToolkit import *

from Bioinformatics_in_Python.utilities import *

with open("Bioinformatics_in_python/test_data/test protein.txt") as f: #opens file
    data = f.read().splitlines()[1:] #reads file, splits lines into list, [1:] removes first line
    dna_seq = "".join(data) #joins list into string without spaces

print(dna_seq)



all_proteins_from_ORF(dna_seq,startReadPos=0, endReadPos=None, ordered=True)


orf_dict = all_proteins_from_ORF(dna_seq, startReadPos=0, endReadPos=None, ordered=True)
for label, orfs in orf_dict.items():
    print(f"{label}:")
    for orf in orfs:
        print(" ", orf)


print(len(dna_seq))
