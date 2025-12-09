import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from bio_seq import bio_seq
from utilities import parse_fasta

# Note: The original file path was hardcoded to a local machine.
# This script now uses a relative path to a test file.
# You may need to change 'test.txt' to the correct FASTA file path.
file_path = "test.txt"

# Transitions and Transversions
# Transition = A <-> G or C <-> T
# Transversion = All other mismatches
seqs = parse_fasta(file_path, seq_type="DNA")

# Ensure you have at least two sequences in the FASTA file
if len(seqs) >= 2:
    seq_1 = seqs[0]
    seq_2 = seqs[1]

    # Use the static method from the bio_seq class
    transitions, transversions, ratio = bio_seq.transitions_transversions_ratio(seq_1, seq_2)

    print(f"Amount of transitions: {transitions}")
    print(f"Amount of transversions: {transversions}")
    if transversions > 0:
        print(f"Ratio (Transitions/Transversions): {ratio:.11f}")
    else:
        print("Ratio is undefined (no transversions).")
else:
    print("Could not find at least two sequences in the provided FASTA file.")
