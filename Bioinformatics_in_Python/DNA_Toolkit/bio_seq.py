from bio_structs import *
from collections import Counter
import random
from utilities import *
import math

class bio_seq:
    """DNA sequnece class. Default value: ATCG, DNA, No label"""

    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label',raise_on_invalid=True):
        """Sequence initialization, validation."""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        if raise_on_invalid == True:
            assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"

    def __validate(self):
        """Check the sequence to make sure it is a valid DNA/RNA/Protein string"""
        return set(BASE[self.seq_type]).issuperset(self.seq) #set command gives a boolean value
    
    def FASTA_input_seq(self, path=None, seq=None, seq_type=None, label=None):
        """Input of a sequence from FASTA format"""
        with open(path, 'r') as f:
            self.label = f.readline().strip()  # Skip the '>' character
            
            remainder = f.read()
            self.seq = remainder.upper().replace("\n", "")

            self.seq_type = input("Enter sequence type (DNA, RNA, Protein): ")
            self.is_valid = self.__validate()
            assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"
            return                

    # @property
    # def seq_len(self):
    #     "return length of sequence"
    #     #need to do this in class, otherwise it won't work outside
    #     return len(self.seq)

    # DNA Toolkit functions:
    def manual_input_seq(self, seq=None, seq_type=None, label=None):
        """Manual input of a sequence"""
        self.seq = input("Enter your sequence: ").upper()
        self.seq_type = input("Enter sequence type (DNA, RNA, Protein): ")
        self.label = input("Enter a label for your sequence: ")
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"


    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type

    def get_seq_len(self):
        """Returns sequence length"""
        return len(self.seq)

    def get_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(BASE[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    def countNucFrequency(self,color=False):
        """"Count frequency of each nucleotide"""
        if color == True:
            return colored(dict(Counter(self.seq)))
        
        else:
            return dict(Counter(self.seq))
         
    # DNA -> RNA
    #make reverse complement of DNA string: DNA to RNA transcription
    def transcription(self):
        """DNA -> RNA Transcription. Replacing Thymine with Uracil"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence"
    
    def reverse_complement(self):
        """Swapping adenine with thymine and guanine with cytosine.
        Reversing newly generated string"""
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC') # maketrans creates a translation table: A->T, C->G, T->A, G->C
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1] # [::-1] reverses the string

    def gc_content_round(self):
        """GC content (rounded) in a DNA/RNA sequence"""
        return f'{round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100)}%'

    #Calculate GC content in %
    def gc_content(self):
        """GC content (.2 decimals) in a DNA/RNA sequence"""
        gc_count = self.seq.count('C') + self.seq.count('G')
        pct = gc_count / len(self.seq) * 100
        return f'{pct:.2f}%'

    def gc_content_subsec(self, k=20, create_file=False):
        """Return a list of GC% strings for each non-overlapping k-mer."""
        results = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i : i + k]
            gc = subseq.count("C") + subseq.count("G")
            pct = gc / len(subseq) * 100
            
            if create_file == True:
                if i == 0:
                    results.append(f"{pct:.2f}% for positions 1-{k}")
                else:
                    results.append(f"{pct:.2f}% for positions {i+1}-{i+k}")
            else:
                results.append(f"{pct:.2f}%")
        return results


    def translate_seq(self, init_pos=0):
        """Translates a DNA sequence into an aminoacid sequence"""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    def codon_usage(self, aminoacid):
        """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        freqDict = dict(Counter(tmpList))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict
    
    def MW_protein_seq(self):
        """Calculates the molecular weight of a protein sequence"""
        if self.seq_type != "Protein":
            print("Not a protein sequence!")

        else:
            return sum(MW_proteins[AA] for AA in self.seq)
        
    def gen_reading_frames(self, startReadPos=0, endReadPos=None):
        """Generates all 6 reading frames of a DNA sequence, including reverse complement"""
        if endReadPos is None or endReadPos > len(self.seq):
            endReadPos = len(self.seq)
        #if endReadPos is not defined or endReadPos is larger than the length of the dna_seq, set endReadPos to the length of the dna_seq
        frames = []
        segment = self.seq[startReadPos:endReadPos] # slice the dna_seq from startReadPos to endReadPos
        dna_seq = segment  # update dna_seq to be the sliced segment
        # Forward frames
        for i in range(3): # 0, 1, 2 positions from which to start translation
            frames.append(self.translate_seq(init_pos=i))
        # Reverse frames
        rev_comp = self.reverse_complement()
        
        original_seq = self.seq  # save the original sequence
        self.seq = rev_comp  # temporarily set self.seq to the reverse complement for translation
        for i in range(3):
            frames.append(self.translate_seq(init_pos=i))
        
        self.seq = original_seq  # restore the original sequence

        return frames

    def all_proteins_from_ORF(self, startReadPos=0, endReadPos=None, ordered=False, create_file = False):
        """
        Return a dict mapping frame labels to lists of every ORF 
        (from 'M' up to but not including the next '_') in that frame. 
        startReadPos and endReadPos represent NUCLEOTIDE positions in the original DNA sequence.
        """
        labels = ['fwd1','fwd2','fwd3','rev1','rev2','rev3']
        startPos = startReadPos
        endPos = endReadPos
        frames = self.gen_reading_frames()
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
        
        if create_file == True:
            file_name = input("Enter file name (e.g. .txt):")
            with open(file_name, "w") as f:
                for label, orfs in orf_dict.items():
                    f.write(f"{label}:\n")
                    for orf in orfs:
                        f.write(f"  {orf}\n")
        
        for label, orfs in orf_dict.items():
            print(f"{label}:")
            for orf in orfs:
                print(" ", orf)
        return orf_dict
    

    def find_motif(self, motif="", amount=0):
        """Find all occurrences of a given motif in a protein sequence and return their starting positions (1-based index).
        In this function, [] only allows for two amino acids, and {} only allows for one amino acid."""
        from bio_structs import BASE
        aa = BASE.get("Protein")
        positions = []
        self.seq = self.seq.upper()
        AA_motif = amount
        motif_check = motif if motif else ""
        consumed = 0
        subsec_motif_correct = True

        if motif == "":
            motif = input("Enter motif (X for any AA, {X} for any AA except X, [XY] for either X or Y): ")
            motif_check = motif
        amount -= 1
        if amount == 0:
            amount = int(input("Enter motif length:")) - 1
            AA_motif = amount + 1

        while len(self.seq) > 0:
            # when motif_check is empty because of a true full match we append; use full_motif_correct flag to be explicit
            if len(motif_check) == 0 and subsec_motif_correct:
                motif_check = motif
                pos = consumed + start + 1
                positions.append(pos)
                self.seq = self.seq[start+1:]
                consumed += start + 1
                continue

            start = self.seq.find(motif[0])
            if start == -1:
                if not positions:
                    print("No motifs found")
                    return
                return positions

            protein_check = self.seq[start:start+AA_motif]

            if motif_check and (motif_check[0] == "X" or motif_check[0] not in aa):
                x = 0

            # assume not full_motif_correct until we prove all tokens full_motif_correct
            full_motif_correct = True
            while len(motif_check) != 0:
                if motif_check[0] == "{":
                    motif_check = motif_check[1:]
                    if not protein_check: #if not tests a condition as true when the list/string is empty. If protein_check would give False (Boolean context)
                        # if we run out of protein to check before motif is consumed, if not function 
                        full_motif_correct = False 
                        motif_check = motif       # mismatch: reset to full motif
                        self.seq = self.seq[start+1:]
                        consumed += start + 1
                        subsec_motif_correct = False
                        break
                    if protein_check[0] != motif_check[0]:
                        motif_check = motif_check[2:]
                        protein_check = protein_check[1:]
                        subsec_motif_correct = True
                        continue
                    else:
                        full_motif_correct = False
                        motif_check = motif       # mismatch: reset
                        self.seq = self.seq[start+1:]
                        consumed += start + 1
                        subsec_motif_correct = False
                        break

                if motif_check[0] == "[":
                    motif_check = motif_check[1:]
                    if not protein_check:
                        full_motif_correct = False
                        motif_check = motif
                        self.seq = self.seq[start+1:]
                        consumed += start + 1
                        subsec_motif_correct = False
                        break
                    if protein_check[0] == motif_check[0] or protein_check[0] == motif_check[1]:
                        motif_check = motif_check[3:]
                        protein_check = protein_check[1:]
                        subsec_motif_correct = True
                        continue
                    else:
                        full_motif_correct = False
                        motif_check = motif
                        self.seq = self.seq[start+1:]
                        consumed += start + 1
                        subsec_motif_correct = False
                        break

                if motif_check[0] == "X":
                    if not protein_check:
                        full_motif_correct = False
                        motif_check = motif
                        self.seq = self.seq[start+1:]
                        consumed += start + 1
                        subsec_motif_correct = False
                        break
                    motif_check = motif_check[1:]
                    protein_check = protein_check[1:]
                    subsec_motif_correct = True
                    continue

                if not protein_check:
                    full_motif_correct = False
                    motif_check = motif
                    self.seq = self.seq[start+1:]
                    consumed += start + 1
                    subsec_motif_correct = False
                    break

                if motif_check[0] == protein_check[0]:
                    motif_check = motif_check[1:]
                    protein_check = protein_check[1:]
                    subsec_motif_correct = True
                    continue

                # literal mismatch
                full_motif_correct = False
                motif_check = motif
                self.seq = self.seq[start+1:]
                consumed += start + 1
                subsec_motif_correct = False
                break

            # if we exited loop and motif_check consumed (len==0) and full_motif_correct True -> record
            if len(motif_check) == 0 and full_motif_correct:
                # set subsec_motif_correct True so the next while-iteration appends and advances
                subsec_motif_correct = True
                continue

        return positions