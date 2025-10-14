
import bio_seq
from bio_structs import *
from utilities import *
import urllib.request

def protein_scrape_uniprot(path="test_data/rosalind_mprt.txt"):
    """Scrape protein sequences from Uniprot given a file with a list of UniProt IDs."""
    with open("test_data/rosalind_mprt.txt", "r") as f:
        temp = [line.strip()[0:6] for line in f]  # Read all lines from the file and strip first 6 characters whitespace/newline characters, you need to call for line in f.readlines() to iterate through each line, else it will treat the whole file as one string
    
    with open("test_data/rosalind_mprt.txt", "r") as f:
        temp2 = [line.strip() for line in f]  # Read all lines from the file and strip first 6 characters whitespace/newline characters, you need to call for line in f.readlines() to iterate through each line, else it will treat the whole file as one string

    temp_dict = {}

    for i, seq in enumerate(temp2):

        var_name = f"{temp2[i]}"  # Create variable names like seq_1, seq_2, etc.
        url_name = f"https://rest.uniprot.org/uniprotkb/{seq[:6]}.fasta"
        temp_dict[var_name] = urllib.request.urlopen(url_name).read().decode('utf-8')

    for i, seq in enumerate(temp_dict.values()):
        
        start = seq.find('\n') + 1     # Update the dictionary to only contain the actual sequence, removing the header line
        temp_dict[f"{temp2[i]}"] = seq[start:]
        temp_dict[f"{temp2[i]}"] = temp_dict[f"{temp2[i]}"].replace("\n","")  
        #cant use var_name here because it will always refer to the last value of i from the previous loop
        #line.replace because strip only removes from start and end of string, not middles of string
    return temp_dict    

proteins = protein_scrape_uniprot()


#N-glycosylation motif: N{P}[ST]{P}
#to allow for the presence of its varying forms, a protein motif is represented by a shorthand as follows: 
# [XY] means "either X or Y" and {X} means "any amino acid except X." 
# For example, the N-glycosylation motif is written as N{P}[ST]{P}.
def find_motif(protein_seq, motif="", amount=0):
    """Find all occurrences of a given motif in a protein sequence and return their starting positions (1-based index).
    In this function, [] only allows for two amino acids, and {} only allows for one amino acid."""
    from bio_structs import BASE
    aa = BASE.get("Protein")
    positions = []
    protein_seq = protein_seq.upper()
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

    while len(protein_seq) > 0:
        # when motif_check is empty because of a true full match we append; use full_motif_correct flag to be explicit
        if len(motif_check) == 0 and subsec_motif_correct:
            motif_check = motif
            pos = consumed + start + 1
            positions.append(pos)
            protein_seq = protein_seq[start+1:]
            consumed += start + 1
            continue

        start = protein_seq.find(motif[0])
        if start == -1:
            if not positions:
                print("No motifs found")
                return
            return positions

        protein_check = protein_seq[start:start+AA_motif]

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
                    protein_seq = protein_seq[start+1:]
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
                    protein_seq = protein_seq[start+1:]
                    consumed += start + 1
                    subsec_motif_correct = False
                    break

            if motif_check[0] == "[":
                motif_check = motif_check[1:]
                if not protein_check:
                    full_motif_correct = False
                    motif_check = motif
                    protein_seq = protein_seq[start+1:]
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
                    protein_seq = protein_seq[start+1:]
                    consumed += start + 1
                    subsec_motif_correct = False
                    break

            if motif_check[0] == "X":
                if not protein_check:
                    full_motif_correct = False
                    motif_check = motif
                    protein_seq = protein_seq[start+1:]
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
                protein_seq = protein_seq[start+1:]
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
            protein_seq = protein_seq[start+1:]
            consumed += start + 1
            subsec_motif_correct = False
            break

        # if we exited loop and motif_check consumed (len==0) and full_motif_correct True -> record
        if len(motif_check) == 0 and full_motif_correct:
            # set subsec_motif_correct True so the next while-iteration appends and advances
            subsec_motif_correct = True
            continue

    return positions

keys = list(proteins.keys())
for i in range(len(keys)):
    key = keys[i]
    seq = proteins[key]
    print(key, end=":\n")
    x = find_motif(seq, motif="N{P}[ST]{P}", amount=4)
    if x:
        print(x) #prevent printing None when no motifs found
    else:
        continue
    write_file("protein_motif_results_test.txt", f">{key}\n{x}", mode="a")  # write results to file

