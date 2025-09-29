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

    def FASTA_input_seq(self, path=None, seq=None, seq_type=None, label=None):
        """Input of a sequence from FASTA format"""
        with open(path, 'r') as f:
            self.label = f.readline().strip()  # Skip the '>' character
            
            remainder = f.read()
            self.seq = remainder.upper().replace("\n", "")

            self.seq_type = input("Enter sequence type (DNA, RNA, Protein): ")
            self.is_valid = self.__validate()
            assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"                

    def parse_fasta(path, seq_type, create_variables=False):
        """
        !! ASSIGN VARIABLE TO FUNCTION TO CREATE LIST !!
        Read a FASTA file with one or more records and returns a list of bio_seq instances, one per record.
        If create_variables is True, create global variables seq_1, seq_2, etc. for each sequence.
        If a sequence is invalid for the given seq_type, it will be skipped with a warning.
        """
         # Avoid circular import
        records = [] #collection of bio_seq objects
        label = None
        seq_parts = [] #list of sequence parts until the next label is found

        with open(path, "r") as f:
            count = 0
            for line in f:
                count += line.count('>')
                
            print(f"There are {count} total sequences in this file.")

        with open(path, "r") as fh:

            
            
            for line in fh: # Iterate through each line in the file
                line = line.strip() # Remove leading/trailing whitespace
                if not line: # Skip empty lines
                    continue

                if line.startswith(">"): #when a new '>' is found, signals start of a new sequence
                    if label is not None: #f label is already set (meaning you just finished reading one record), you:
                        full_seq = "".join(seq_parts) #join all parts of the sequence into one string
                        temp_seq = bio_seq(full_seq, seq_type, label,raise_on_invalid=False) #create a new bio_seq object and add it to the records list
                        if temp_seq.is_valid == True:
                            records.append(temp_seq) # add the bio_seq object to the records list only if it's valid
                        else:
                            print(f"Warning: Skipping invalid {seq_type} sequence with label '{label}'")
                    label = line[1:] #if label is None, set label to the current line without '>'
                    seq_parts = [] #reset seq_parts for the next sequence
                else:
                    seq_parts.append(line.upper()) # add the current line (part of the sequence) to seq_parts, converting to uppercase

        if label is not None: #read the last record after the loop ends
            full_seq = "".join(seq_parts)
            temp_seq = bio_seq(full_seq, seq_type, label,raise_on_invalid=False)
            if temp_seq.is_valid:
                records.append(temp_seq)
            else:
                print(f"Warning: Sequence '{label}' is invalid for type '{seq_type}' and will be skipped.")
        print(f"Successfully parsed {len(records)} valid {seq_type} sequences from the file.")
        
        if create_variables == True:

            #The enumerate function wraps any iterable and yields pairs of (index, element), allowing you to track positions without manual counter management.
            for i, record in enumerate(records): # Loop through each record with its index, index starts at 0 and ends at records length -1
                var_name = f"seq_{i+1}" # Create a variable name like seq_1, seq_2, etc.
                globals()[var_name] = record #globals() returns the dictionary representing the current module’s global variables.
                                            #Assigning to globals()[var_name] creates a new global variable with the name in var_name.
                print(f"Created variable '{var_name}' for sequence with label '{record.label}'")

        else:
            return records
    

    # DNA Toolkit functions:
    def manual_input_seq(self, seq=None, seq_type=None, label=None):
        """Manual input of a sequence"""
        self.seq = input("Enter your sequence: ").upper()
        self.seq_type = input("Enter sequence type (DNA, RNA, Protein): ")
        self.label = input("Enter a label for your sequence: ")
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"

    def __validate(self):
        """Check the sequence to make sure it is a valid DNA/RNA/Protein string"""
        return set(BASE[self.seq_type]).issuperset(self.seq) #set command gives a boolean value

    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type

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

    def gc_content_subsec(self, k=20):
        """Return a list of GC% strings for each non-overlapping k-mer."""
        results = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i : i + k]
            gc = subseq.count("C") + subseq.count("G")
            pct = gc / len(subseq) * 100
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

    def all_proteins_from_ORF(self, startReadPos=0, endReadPos=None, ordered=False):
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
        return orf_dict