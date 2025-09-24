from Bioinformatics_in_Python.DNAToolkit import *

from Bioinformatics_in_Python.DNAToolkit import translate_seq, codon_usage

codon = translate_seq(dna_seq, init_pos=0)
codon = ''.join(codon)  # Join the list of codons into a string
codon

from Bioinformatics_in_Python.utilities import read_file, FASTA_to_dict, gc_content, colored