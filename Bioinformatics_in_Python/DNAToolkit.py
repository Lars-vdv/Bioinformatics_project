import random
from Bioinformatics_in_Python.structures import *
from Bioinformatics_in_Python.utilities import *


def generate_dna_sequence(length):
    nuc = 'ATCG'
    return ''.join(random.choice(nuc) for _ in range(length))

dna_seq = generate_dna_sequence(100)  # Generate a random DNA sequence of length 100


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
    return "  ".join(f"{nuc}:{tmpFreqDict[nuc]}" for nuc in "ACGT")

print(f"{colored(countNucFrequency(dna_seq))}")#in color, remove color for normal
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
    return ''.join([reverse_dna_seq[nuc] for nuc in dna_seq])  # swaps each nucleotide for its complement

print(complement(dna_seq))

# Reverse complement of DNA string
def reverse_complement(dna_seq):
    """Creates reverse complement of default DNA string"""
    return ''.join([reverse_dna_seq[nuc] for nuc in dna_seq][::-1])  # swaps and then reverses

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

#Calculate rounded GC content in %
def gc_content_round(dna_seq):
    """GC content (rounded) in a DNA/RNA sequence"""
    return f'{round((dna_seq.count('C') + dna_seq.count('G')) / len(dna_seq) * 100)}%'
                 #round -> rounds number, count -> counts amount of how much a defined variable is in the string, len -> counts lenght of new string (G+C) *100 makes percentage
                #function has to be inside {} for f string. characters (like%) outside. All has to between f'..'

#gc_content(dna_seq)

#Calculate GC content in %
def gc_content(dna_seq):
    """GC content (.2 decimals) in a DNA/RNA sequence"""
    gc_count = dna_seq.count('C') + dna_seq.count('G')
    pct = gc_count / len(dna_seq) * 100
    return f'{pct:.2f}%'
                 #round -> rounds number, count -> counts amount of how much a defined variable is in the string, len -> counts lenght of new string (G+C) *100 makes percentage
                #function has to be inside {} for f string. characters (like%) outside. All has to between f'..'
                #:.2f - > : starts format specification, .2 means two digits, f means fixed-point notation

#Return a list of GC% strings for each non-overlapping k-mer
def gc_content_subsec(dna_seq, k=20):
    """
    Return a list of GC% strings for each non-overlapping k-mer.
    """
    #return [....] creates list
    #[expression for variable in iterable]
    return f'{round(
        gc_content(dna_seq[i : i + k]) #dna_seq[start : end], When there are no more i values in the range, the loop ends.
        for i in range(0, len(dna_seq) - k + 1, k) #range(start, stop, step)
        # range start=0, stop=len(dna_seq)-k+1 ensures full windows only, step=k, +1 is needed because python excludes stop index 
        #len(dna_seq) - k is the highest starting index that still leaves k-amount characters to slice.
    )}%'



#Create dictonary for FASTA file
def FASTA_to_dict(file_path):
    """
    Convert a FASTA file (path) into a dictionary.
    """
    FASTAfile = read_file(file_path)
    FASTAdict = {}
    FASTAlabel = ""
    
    for line in FASTAfile:
        if line.startswith('>'):
            FASTAlabel = line.strip()  # Remove newline and spaces
            FASTAdict[FASTAlabel] = ""  # Initialize the key in the dictionary
        else:
            FASTAdict[FASTAlabel] += line.strip()  # Append sequence to the current label
    
    return (
        FASTAdict
    )



#Create GC content dictionary
def GC_dict(FASTAdict,ndigits=2):
    """Returns a dictionary with rounded GC content for each sequence in FASTAdict"""
    temp = {
        key: round(gc_content(value),ndigits) for (key, value) in FASTAdict.items()
        } #substitute DNA sequence from FASTAdict with GC content, FASTAdict.items() returns a list of tuples (key, value) for each item in the dictionary, key: FASTAlabel, value: DNA sequence
    return temp

#Create GC content_subseq dictionary
def GC_dict_subseq(FASTAdict,ndigits=2):
    """Returns a dictionary with GC content_subseq for each sequence in FASTAdict"""
    temp = {key: round(gc_content_subseq(value),ndigits) for (key, value) in FASTAdict.items()} #substitute DNA sequence from FASTAdict with GC content, FASTAdict.items() returns a list of tuples (key, value) for each item in the dictionary, key: FASTAlabel, value: DNA sequence
    return temp

#Find max and min GC content
#maxGCkey = max(GC_dict, key=GC_dict.get) #max() returns the key with the highest value in GC_dict, key=GC_dict.get tells max() to use the values of GC_dict to compare
#print(f"The sequence with the highest GC content is {maxGCkey} with a GC content of {GC_dict[maxGCkey]}%") #f-string to format the output
#maxGCkey and minGCkey .
def max_min_GCkey(GC_dict):
    """Returns the keys with the highest and lowest GC content"""
    maxGCkey = max(GC_dict, key=GC_dict.get)
    minGCkey = min(GC_dict, key=GC_dict.get)
    print(f"The sequence with the highest GC content is {maxGCkey} with a GC content of {GC_dict[maxGCkey]}%")
    print(f"The sequence with the lowest GC content is {minGCkey} with a GC content of {GC_dict[minGCkey]}%")

#Translate DNA sequence to amino acid sequence
def translate_seq(dna_seq, init_pos=0):
    """Translates a DNA sequence into an amino acid sequence starting at init_pos"""
    return (
        [DNA_Codons[dna_seq[pos:pos + 3]] for pos in range(init_pos, len(dna_seq) - 2, 3)]
    )

#frequency of each codon encoding a given aminoacid in a DNA sequence
def codon_usage(dna_seq, aminoacid):
    """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
    from collections import Counter
    tmpList = []
    for i in range(0, len(dna_seq) - 2, 3):
        if DNA_Codons[dna_seq[i:i + 3]] == aminoacid:
            tmpList.append(dna_seq[i:i + 3])

    freqDict = dict(Counter(tmpList))
    totalWight = sum(freqDict.values())
    for dna_seq in freqDict:
        freqDict[dna_seq] = round(freqDict[dna_seq] / totalWight, 2)
    return freqDict