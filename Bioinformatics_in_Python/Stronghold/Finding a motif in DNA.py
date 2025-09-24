from Bioinformatics_in_Python.DNAToolkit import *

import re

with open("Bioinformatics_in_Python/test_data/rosalind_subs.txt") as f:
    DNA_seq = f.readline().strip()   # reads line 1
    motif   = f.readline().strip()   # reads line 2


def find_motif(dna_seq, motif):
    """Find all occurrences of a motif in a DNA sequence."""
    pattern = f'(?={motif})'  # Lookahead assertion to find overlapping motifs
    positions = [m.start() + 1 for m in re.finditer(pattern, dna_seq)]  # +1 to convert to 1-based index
    return positions

x = find_motif(DNA_seq, motif)
print(" ".join(map(str, x)))  # Print positions as a space-separated string
# The map function converts each integer in the list to a string for joining
# The join function combines the strings with a space in between