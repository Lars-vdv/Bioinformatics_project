import os
print(os.getcwd())
print(os.path.exists("test_data/rosalind_lcsm.txt"))

from utilities import *
from bio_structs import *
from bio_seq import *
from difflib import SequenceMatcher

seq = parse_fasta(path=r"test_data\rosalind_lcsm.txt", seq_type="DNA", create_variables=False)

for i in range(len(seq)):
    print(getattr(seq[i], "seq_len"))

print(getattr(seq[0], "seq_len")) #test if parsing was successful

# for i in len(getattr(seq, "seq")):
#     seq_strings = []
#     seq__strings = seq[i]

seq = sorted(seq, key=lambda s: s.seq_len, reverse=True)

# #first order the sequences on length
# seq = sorted(seq, key=lambda s: , reverse=True)
# print(seq)

for i in len(seq):
    print(getattr(seq[i], "seq_len"))



common_substrings = {}
for i in range(len(seq)):
    for j in range(i+1, len(seq)):
        string1 = getattr(seq[i], "seq")
        string2 = getattr(seq[j], "seq")
        match = SequenceMatcher(None, string1, string2).find_longest_match(0, len(string1), 0, len(string2)) 
        #find.longest_match is based on a[start1:end1] == b[start2:end2]
        matching_substring = string1[match.a: match.a + match.size] # Extract the matching substring
        if matching_substring not in common_substrings:
            common_substrings[matching_substring] = 1
        else:
            common_substrings[matching_substring] += 1

print(common_substrings)

max_len = max(len(k) for k in common_substrings) 
#max finds the longest length of the keys in the dictionary
#len(k) for k in common_substrings iterates over all the keys of the dictionary. You dont need a forloop nor keys() method
lcs = [k for k in common_substrings if len(k) == max_len] 
print(lcs)
#list comprehension to create a list of all keys with the max length
#k — the value yielded for the resulting list (each matched element).
#for k in common_substrings — iterate over every item in the iterable common_substrings and bind it to k.
#if len(k) == max_len — include k in the output only when this boolean expression is true.