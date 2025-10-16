import os
from utilities import *
from bio_structs import *
from bio_seq import *


# def longest_common_substring(fasta_file_path):
#     """
#     Finds the longest common substring in a FASTA file.
#     """
#     seqs = parse_fasta(path=fasta_file_path, seq_type="DNA", create_variables=False)

#     if not seqs:
#         return ""

#     # Sort sequences by length and get the shortest one
#     seqs = sorted(seqs, key=lambda s: s.seq_len)
#     shortest_seq = seqs[0].seq
#     other_seqs = [s.seq for s in seqs[1:]]

#     # Iterate through all possible substrings of the shortest sequence, from longest to shortest
#     for length in range(len(shortest_seq), 0, -1):
#         for i in range(len(shortest_seq) - length + 1):
#             substring = shortest_seq[i:i+length]
            
#             # Check if this substring exists in all other sequences using the optimized `in` operator
#             if all(substring in other for other in other_seqs):
#                 return substring
                           

# Find the longest common substring
lcs = longest_common_substring(fasta_file_path=r"test_data\rosalind_lcsm.txt")
