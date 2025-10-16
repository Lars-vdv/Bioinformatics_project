import os
from utilities import *
from bio_structs import *
from bio_seq import *
from difflib import SequenceMatcher

def longest_common_substring(fasta_file):
    """
    Finds the longest common substring in a FASTA file.
    """
    seqs = parse_fasta(path=fasta_file, seq_type="DNA", create_variables=False)

    if not seqs:
        return ""

    # Sort sequences by length and get the shortest one
    seqs = sorted(seqs, key=lambda s: s.seq_len)
    shortest_seq = seqs[0].seq
    other_seqs = [s.seq for s in seqs[1:]]

    # Iterate through all possible substrings of the shortest sequence, from longest to shortest
    for length in range(len(shortest_seq), 0, -1):
        for i in range(len(shortest_seq) - length + 1):
            substring = shortest_seq[i:i+length]

            # Check if this substring exists in all other sequences
            if all(substring in other for other in other_seqs):
                return substring

    return ""

if __name__ == "__main__":
    # Construct the full path to the file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, "test_data", "rosalind_lcsm.txt")

    # Find the longest common substring
    lcs = longest_common_substring(file_path)

    # Print the result
    print(f"Longest common substring: {lcs}")