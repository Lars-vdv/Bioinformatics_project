import os
print(os.getcwd())


from Bioinformatics_in_Python.DNAToolkit import *

from Bioinformatics_in_Python.utilities import *

with open("Bioinformatics_in_Python/test_data/rosalind_hamm.txt") as f:
    dna_seq1 = f.readline().strip()  # Read the first line and strip any whitespace/newline characters
    dna_seq2 = f.readline().strip()  # Read the second line and strip any whitespace/newline characters

print(dna_seq1)
print(dna_seq2)

#to create a dictionary instead of two separate variables:
    # dna_seqs = []
    # with open("Bioinformatics_in_python/test_data/test protein.txt") as f:
    #     for line in f:
    #         dna_seqs.append(line.strip('\n'))

# positions where characters differ
mismatch_positions = [
    i for i, (a, b) in enumerate(zip(dna_seq1, dna_seq2)) if a != b #enumerate gives index and value, zip pairs elements from both sequences
] 
# i for i, (a,b) means “unpack each item from the iterable into two variables, then yield the variable i.”

# number of mismatches
mismatch_count = len(mismatch_positions)

print(mismatch_positions)  # e.g. [3, 4]
print(mismatch_count)      # e.g. 2

mismatches = [
    (i, a, b) for i, (a, b) in enumerate(zip(dna_seq1, dna_seq2)) if a != b
]
print(mismatches)

