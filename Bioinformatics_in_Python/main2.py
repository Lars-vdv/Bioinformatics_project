from bio_seq import *
from bio_structs import *
from utilities import *
from pprint import pprint
import math

sequences_data = parse_fasta(path="test_data/rosalind_gc.txt", seq_type="DNA", create_variables=True)

print(sequences_data['seq_1'].get_seq_info()) #alternative
globals().update(sequences_data) #see OneNote for explanation
print(seq_1.get_seq_info())

seq_1dict = seq_1.all_proteins_from_ORF(ordered=True)


write_file("test_output.txt", seq_1.seq)

for rf in seq_1.all_proteins_from_ORF():

    write_file("test_output.txt", rf)
# for rf in seq_1.gen_reading_frames():
#     write_file("test_output.txt", rf)

#Check if the returned data is a dictionary and not empty
# if isinstance(sequences_data, dict) and sequences_data:
#     # This is the correct way to dynamically create global variables in the main script's scope.
#     globals().update(sequences_data)
    
#     # Now, variables like seq_1, seq_2, etc., exist and can be accessed directly.
#     # The script will fail if 'seq_1' was not in the parsed data, so we check for its existence.
#     if 'seq_1' in globals():
#         print(seq_1.get_seq_info())
#         print(seq_1.countNucFrequency(color=True))
#         print(colored(seq_1.transcription()))
#         print(seq_1.reverse_complement())
#         print(seq_1.gc_content_round())
#         print(seq_1.gc_content())
#         print(seq_1.gc_content_subsec())
#         print(seq_1.translate_seq())
#         print(seq_1.codon_usage('L'))
#         print(seq_1.gen_reading_frames())
#         print(seq_1.all_proteins_from_ORF(ordered=True))
#     else:
#         print("Variable 'seq_1' was not created. Check the input file.")
# else:
#     print("No valid sequences were found or the data format was incorrect.")
