from bio_seq import *
from bio_structs import *
from utilities import *
from pprint import pprint
import math

#set variable (e.g. random_sequence) as an instance of bio_seq class
#generate random DNA sequence of length 100
random_sequence =bio_seq()
random_sequence.generate_rnd_seq(length=100, seq_type="DNA")


########## general sequence operations ##########
#get sequence information
print(random_sequence.get_seq_info())

# count nucleotide frequency
print(random_sequence.countNucFrequency(color=True))

#transcribe DNA sequence into (colored) RNA sequence
print(colored(random_sequence.transcription()))

#DNA reverse complement,
print(random_sequence.reverse_complement())     



# GC content (rounded), 
print(random_sequence.gc_content_round())

#GC content (float),
print(random_sequence.gc_content())

#GC content in subsections (k=10),
print(random_sequence.gc_content_subsec())

#  translate DNA/RNA sequence,
print(random_sequence.translate_seq())

#  codon usage for e.g. Leucine,
print(random_sequence.codon_usage('L'))

# generate  all reading frames from sequence,
print(random_sequence.gen_reading_frames())

#generate all proteins from ORF (ordered)
print(random_sequence.all_proteins_from_ORF(ordered=True))

#find motif in DNA/RNA sequence (e.g. ATG)
print(random_sequence.find_DNA_RNA_motif(motif="ATG"))


########## Protein specific operations ##########
#generate random Protein sequence of length 50
random_protein_sequence = bio_seq()
random_protein_sequence.generate_rnd_seq(length=50, seq_type="Protein")

#generate molecular weight of protein
print(MW_proteins(random_protein_sequence.seq))

#find motif in protein sequence (e.g. N{P}[ST]{P})
print(random_protein_sequence.find_motif(motif="N{P}[ST]{P}"))

########## static methods (no bio_seq object needed) ##########

#find longest common substring in fasta file e.g. test_files/longest_common_substring.fasta
#use bio_seq.xxx to call static method
fasta_file_path = "test_files/longest_common_substring.fasta"
print(bio_seq.longest_common_substring(fasta_file_path))

#splice gene and translate to protein from fasta file e.g. test_files/splice.fasta
#or use dictionary with gene and introns
#genes_and_introns should be a dictionary with complete gene[0] and intron sequences[1:n]

file_path = "test_files/splice.fasta"
gene_and_introns = bio_seq(seq="ACTG", seq_type="DNA")
gene_and_introns = parse_fasta(path = file_path, seq_type="DNA", create_variables=True)
v = list(gene_and_introns.values())

print(gene_and_introns["seq_1"]) #print whole gene (introns + exons))
#get the sequence of the objects in the dictionary
for i in range(len(gene_and_introns)):
    name = str(f"seq_{i+1}")
    
    gene_and_introns[name] = getattr(gene_and_introns[name],"seq")

bio_seq.splice_gene(gene_and_introns)


