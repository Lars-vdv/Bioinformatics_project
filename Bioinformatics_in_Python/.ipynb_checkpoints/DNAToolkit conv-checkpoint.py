#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Converted from Jupyter Notebook: notebook.ipynb
Conversion Date: 2025-08-05T13:13:48.455Z
"""

from utilities import *

from IPython.display import display, HTML

from structures import * # * -> import all


dna_seq = "ATCGATCGTTACCAGGGCCCCcgatcgatcgattcagcaatacccggggcatacgactagctagcgagCTAGCGTAGCATGACTGACTGATCG"

# DNA toolset/code testing file
def validateSeq(dna_seq): # def defines a function, in this case validateSeq
    """"Validate if sequence is DNA sequence"""
    tmpseq = dna_seq.upper() # create temporary sequence file, .upper makes everything uppercase
    for nuc in tmpseq:  #nuc is created here, it loops througyh every character and names it nuc, which you can use as 'x'
        if nuc not in nucleotides:
            return False
    return tmpseq

dna_seq = validateSeq(dna_seq)
print(dna_seq)

# Count nucleotide frequency
def countNucFrequency(dna_seq):
    """"Count frequency of each nucleotide"""
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0} #{} create temporary dictonary
    for nuc in dna_seq:
        tmpFreqDict[nuc] += 1 #For every nucleotide counted add +1, nuc is added to dictonary value because it is the only value in the dictonary
    return tmpFreqDict

print(f'{countNucFrequency(dna_seq)}') #in color, remove color for normal
###print(*map(dna_seq.count, "ACGT")) #also works

# DNA -> RNA
#make reverse complement of DNA string: DNA to RNA transcription
def transcription(dna_seq):
    """"Transcribes default DNA string to RNA"""
    return dna_seq.replace("T","U") #here the string is returned in RNA seq, .replace replaces the T's with U's 

rna_seq = transcription(dna_seq)
print(rna_seq)

# NOT REVERSE Complement for DNA string 
def complement(dna_seq):
    """"Creates complement of default DNA string"""
    return ''.join([reverse_dna_seq[nuc] for nuc in dna_seq]) # [reverse..][nuc] swaps the letter for every nuc in dna_seq (so every letter) as defined in reverse_dna_seq dict
    
print(complement(dna_seq))

# Reverse complement for DNA string 
def reverse_complement(dna_seq):
    """"Creates reverse complement of default DNA string"""
    return ''.join([reverse_dna_seq[nuc] for nuc in dna_seq])[::-1] #[::-1] reverses the string, [reverse..][nuc] swaps the letter for every nuc in dna_seq (so every letter) as defined in reverse_dna_seq dict
    
print(reverse_complement(dna_seq))

#Print all sequences and transcript
def print_all_seq(dna_seq):
    """"Prints all DNA and RNA string versions"""
    print(f"DNA sequence: {dna_seq}\n\n"
      f"Complement DNA sequence: {complement(dna_seq)}\n\n"
      f"Reverse complement DNA sequence: {reverse_complement(dna_seq)}\n\n"
      f"RNA transcript of DNA: {rna_seq[::-1]}\n\n"
      f"RNA transcript of complement DNA: {rna_seq}\n")

print_all_seq(dna_seq)

#Print all sequences and transcript in colors
def print_all_seq_colored(dna_seq):
    """"Prints all DNA and RNA string versions in colors"""
    print(f"DNA sequence: {colored(dna_seq)}\n\n"
      f"Complement DNA sequence: {colored(complement(dna_seq))}\n\n"
      f"Reverse complement DNA sequence: {colored(reverse_complement(dna_seq))}\n\n"
      f"RNA transcript of DNA: {colored(rna_seq[::-1])}\n\n"
      f"RNA transcript of complement DNA: {colored(rna_seq)}\n\n")

print_all_seq_colored(dna_seq)

print(dna_seq)

#Calculate GC content in %
def gc_content(seq):
    """GC content in a DNA/RNA sequence"""
    return f'{round((seq.count('C') + seq.count('G')) / len(seq) * 100)}%'
                 #round -> rounds number, count -> counts amount of how much a defined variable is in the string, len -> counts lenght of new string (G+C) *100 makes percentage
                #function has to be inside {} for f string. characters (like%) outside. All has to between f'..'

#gc_content(dna_seq)

def gc_content_subsec(seq, k=20):
    """
    Return a list of GC% strings for each non-overlapping k-mer.
    """
    #return [....] creates list
    #[expression for variable in iterable]
    return [
        gc_content(seq[i : i + k]) #seq[start : end], When there are no more i values in the range, the loop ends.
        for i in range(0, len(seq) - k + 1, k) #range(start, stop, step)
        # range start=0, stop=len(seq)-k+1 ensures full windows only, step=k, +1 is needed because python excludes stop index 
        #len(seq) - k is the highest starting index that still leaves k-amount characters to slice.
    ]

gc_content_subsec(dna_seq,k=10)

pwd()