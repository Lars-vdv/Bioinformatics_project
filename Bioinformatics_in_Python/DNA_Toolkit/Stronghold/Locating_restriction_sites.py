import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from bio_seq import bio_seq
from bio_structs import *
from utilities import *

revp = bio_seq()
revp = parse_fasta(r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\test_data\rosalind_revp.txt", seq_type="DNA", create_variables=True)

seq_1 = revp.get("seq_1")
print(seq_1)


seq = getattr(seq_1, "seq")
print(seq)

def find_restriction_sites(DNA_seq, min_length=4, max_length=12):
    """This function can find possible restriction sites based on reverse palindromes
    of a DNA string"""
    DNA_seq = bio_seq(DNA_seq)
    comp_DNA_seq = DNA_seq.reverse_complement()
    length = getattr(DNA_seq, "seq_len")
    DNA_seq = getattr(DNA_seq, "seq")
    dict = {"location":"length"}
    for i in range(length):
        for j in range(min_length, max_length+1):
            if i + j <= length:
               fragment = DNA_seq[i:i+j]
               rev_comp_fragment = comp_DNA_seq[length - (i+j):length - i]
               if fragment == rev_comp_fragment:
                    dict.update({i : j})

    return dict


test = find_restriction_sites(DNA_seq=seq)

write_file("test.txt", test)