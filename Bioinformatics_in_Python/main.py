from Bioinformatics_in_Python.DNAToolkit import *

with open("Bioinformatics_in_python/test_data/test protein.txt") as f: #opens file
    data = f.read().splitlines()[1:] #reads file, splits lines into list, [1:] removes first line
    dna_seq = "".join(data) #joins list into string without spaces

# print(f'\nSequence: {dna_seq}\n')

# print(f'[1] + Sequence Length: {len(dna_seq)}\n')

# print(f'[2] + Nucleotide Frequency: {countNucFrequency(dna_seq)}\n')

# print(f'[3] + DNA/RNA Transcription: {transcription(dna_seq)}\n')

# print(f"[4] + DNA String + Complement + Reverse Complement:\n5' {dna_seq} 3'")

# print(f"   {''.join(['|' for c in range(len(dna_seq))])}")

# print(f"3' {complement(dna_seq)} 5' [Complement]")

# print(f"5' {reverse_complement(dna_seq)} 3' [Rev. Complement]\n")

# print(f'[5] + GC Content: {gc_content(dna_seq)}\n')

# print(f'[6] + GC Content in Subsection k=5: {gc_content_subseq(dna_seq, k=5)}\n')

# print(f'[7] + Aminoacids Sequence from DNA: {translate_seq(dna_seq, 0)}\n')


# print(f'[8] + Codon frequency (L): {codon_usage(dna_seq, "L")}\n')


# print('[9] + Reading_frames:')
# for frame in gen_reading_frames(dna_seq):
#     print(frame)

def all_proteins_from_ORF(dna_seq, startReadPos=0, endReadPos=None, ordered=False):
    """
    Return a dict mapping frame labels to lists of every ORF 
    (from 'M' up to but not including the next '_') in that frame. 
    startReadPos and endReadPos represent NUCLEOTIDE positions in the original DNA sequence.
    """
    labels = ['fwd1','fwd2','fwd3','rev1','rev2','rev3']
    startPos = startReadPos
    endPos = endReadPos
    frames = gen_reading_frames(dna_seq, startReadPos=startPos, endReadPos=endPos)
    orf_dict = {}

    for label, frame in zip(labels, frames): # loops through each frame with its corresponding label
        prot = ''.join(frame)
        orfs = []
        search_pos = 0

        while True: # loop to find all ORFs in the protein sequence
            # find the next 'M'
            start = prot.find('M', search_pos)
            if start == -1:
                break

            # find the next '_' after that 'M'
            stop = prot.find('_', start + 1)
            if stop == -1:
                break

            orfs.append(prot[start:stop])
            # advance past this 'M' so we can catch overlapping ORFs 
            search_pos = stop + 1 # move search position to just after the found stop codon to prevent duplicates

        orf_dict[label] = orfs

    if ordered:
        for label in orf_dict: #here you loop through the keys of the dictionary = label, otherwise only last label would be sorted
            orf_dict[label].sort(key=len, reverse=True) # Sort ORFs by length in descending order if ordered is True
    return orf_dict
    
orf_dict = all_proteins_from_ORF(dna_seq, startReadPos=0, endReadPos=None, ordered=True)
for label, orfs in orf_dict.items():
    print(f"{label}:")
    for orf in orfs:
        print(" ", orf)