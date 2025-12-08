import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

print(os.getcwd())
file_path = r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\test_data\rosalind_splc.txt"

from bio_seq import bio_seq
from bio_structs import *
from utilities import *
import re

gene_and_introns = bio_seq(seq="ACTG", seq_type="DNA")
gene_and_introns = parse_fasta(path = file_path, seq_type="DNA", create_variables=True)


v = list(gene_and_introns.values())

print(gene_and_introns["seq_1"])

#get the sequence of the objects in the dictionary
for i in range(len(gene_and_introns)):
    name = str(f"seq_{i+1}")
    
    gene_and_introns[name] = getattr(gene_and_introns[name],"seq")



def find_motif(dna_seq, motif):
    """Find all occurrences of a motif in a DNA sequence."""
    pattern = f'(?={motif})'  # Lookahead assertion to find overlapping motifs
    position_start = [m.start() + 1 for m in re.finditer(pattern, dna_seq)]  # +1 to convert to 1-based index
    return position_start

for j in range(len(gene_and_introns)):
    
    #set whole gene sequence
    if j == 0:
        gene = gene_and_introns.get("seq_1")

    #once first intron is removed, use new gene generated from removing the intron
    #repeat this
    else:
        gene = gene_without_introns

    #set name of key in dictionary so you can retrieve it
    val = f"seq_{j+2}"

    #get value of key
    intron = gene_and_introns.get(val)

    #last round no intron, break
    if intron == None:
        break

    length_intron = len(intron)

    #find the intron location in gene, start position is +1 indexing (See find motif function)
    position_start = find_motif(gene, intron)
    position_start = position_start[0]

    #position_end is the start  position + the length of the intron
    position_end = position_start + length_intron

    #slice out the intron
    gene_without_introns = gene[:position_start] + gene[position_end:]


# print(gene_without_introns)
# print(len(gene_without_introns))
# print(gene_and_introns["seq_1"])
# print(len(gene_and_introns["seq_1"]))

#convert sequence to bio.seq object and translate it into protein sequence
protein = bio_seq(seq=gene_without_introns)
protein = bio_seq.translate_seq(protein)
protein = "".join(protein)
print(protein)



write_file("splice.txt", protein)

bio_seq.splice_gene(gene_and_introns)