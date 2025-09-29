from bio_seq import bio_seq
from bio_structs import *
from utilities import *
from pprint import pprint
import math

               
parse_fasta(path=r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\test_data\rosalind_gc.txt", seq_type="DNA", create_variables=True)


print(seq_1.get_seq_info())
print(seq_1.countNucFrequency(color=True))
print(colored(seq_1.transcription()))
print(seq_1.reverse_complement())
print(seq_1.gc_content_round())
print(seq_1.gc_content())
print(seq_1.gc_content_subsec())
print(seq_1.translate_seq())
print(seq_1.codon_usage('L'))
print(seq_1.gen_reading_frames())
print(seq_1.all_proteins_from_ORF(ordered=True))


