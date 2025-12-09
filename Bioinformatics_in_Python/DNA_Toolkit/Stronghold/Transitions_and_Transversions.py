import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

print(os.getcwd())
file_path = r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\test_data\rosalind_tran.txt"

from bio_seq import bio_seq
from bio_structs import *
from utilities import *

#Transitions and Transversions
#Transisiton = A <-> G or C <-> T
#Transversion = A > C/T, G > C/T, T -> A/G 
seqs = parse_fasta(file_path, seq_type="DNA", create_variables=True)

seq_1 = seqs.get("seq_1")
seq_2 = seqs.get("seq_2")

# def transitions_transversions_ratio(seq_1: bio_seq, seq_2: bio_seq):
#     Transitions = 0
#     Transversions = 0

#     for a, b in zip(seq_1.seq, seq_2.seq):
#         if a != b:
#             if (a == "A" and b == "G") or (a == "G" and b == "A") or (a == "C" and b == "T") or (a == "T" and b == "C"):
#                 Transitions += 1
#             else:
#                 Transversions += 1
#         else:
#             continue
#     print(f"Amount of transitions: {Transitions}\nAmount of transversions: {Transversions}\nRatio (Transitions/Transversions): {Transitions/Transversions}")


def transitions_transversions_ratio(seq_1, seq_2):
    Transitions = 0
    Transversions = 0

    bool_1 = isinstance(seq_1, bio_seq)
    if bool_1 == True:
        seq_1 = getattr(seq_1, 'seq')
    
    bool_2 = isinstance(seq_2, bio_seq)
    if bool_2 == True:
        seq_2 = getattr(seq_2, 'seq')

    for a, b in zip(seq_1, seq_2):
        if a != b:
            if (a == "A" and b == "G") or (a == "G" and b == "A") or (a == "C" and b == "T") or (a == "T" and b == "C"):
                Transitions += 1
            else:
                Transversions += 1
        else:
            continue

    if Transversions == 0:
        print("No transversions found, ratio is undefined (division by zero).")
        return
    
    if Transitions == 0:
        print("No transitions found, ratio is 0.")
        return

    print(f"Amount of transitions: {Transitions}\nAmount of transversions: {Transversions}\nRatio (Transitions/Transversions): {Transitions/Transversions}")

transitions_transversions_ratio(seq_1, seq_2)

