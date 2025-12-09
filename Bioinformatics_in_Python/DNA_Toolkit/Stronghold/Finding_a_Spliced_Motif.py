import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

print(os.getcwd())
file_path = r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\test_data\rosalind_sseq.txt"

from bio_seq import bio_seq
from bio_structs import *
from utilities import *

#parse fasta and assign sequences
seqs = parse_fasta(file_path, "DNA")
seq = seqs[0]
motif = seqs[1]

print(getattr(motif, 'seq'))

#find spliced motif
#rosalind function. bio_seq only.
#FIND BIO_SEQ AND RAW STRing function below and in bio_seq.py
def find_spliced_motif_rosalind(seq: bio_seq, motif: bio_seq):
    """Finds the positions of a spliced motif in a sequence."""
    """bio_seq only function (rosalind problem)"""
    index = []
    motif_length = len(motif.seq)
    motif_position = 0
    for i in range(len(seq.seq)):
        if seq.seq[i] == motif.seq[motif_position]:
            index.append(i+1)
            motif_position += 1

        else:
            continue
        
        if motif_position == motif_length:
            print("Positions found: ", index)
            return index

pos = find_spliced_motif_rosalind(seq, motif)

seq1 = bio_seq()
seq2 = bio_seq()


seq1.generate_rnd_seq(200, "DNA")
seq2.generate_rnd_seq(100, "DNA")

print(bio_seq.find_spliced_motif(seq1, seq2))


# list_as_string = ' '.join(map(str, find_spliced_motif(seq, motif)))
# print(list_as_string)
