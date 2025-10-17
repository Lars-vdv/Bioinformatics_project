import sys
import urllib.request

def colored(data):
    """Return a colored string for a sequence or a dict of nucleotide counts."""
    colors = {
        'A': '\033[31m',  # red
        'C': '\033[34m',  # blue
        'G': '\033[33m',  # yellow
        'T': '\033[32m',  # green
        'U': '\033[35m',  # magenta
        'reset': '\033[0m'
    }

    # If it's a dict, build a "{ A:1, C:2, ... }" string with colored keys
    if isinstance(data, dict):
        items = []
        for nuc, count in data.items():
            code = colors.get(nuc, colors['reset'])
            items.append(f"{code}{nuc}{colors['reset']}: {count}")
        return "{ " + ", ".join(items) + " }"

    # Otherwise assume it's a sequence string
    result = ""
    for nuc in data:
        code = colors.get(nuc, colors['reset'])
        result += f"{code}{nuc}{colors['reset']}"
    return result


def colored_html(seq):
    colors = {
        'A': 'red',
        'C': 'blue',
        'G': 'orange',
        'T': 'green',
        'U': 'magenta'
    }
    html = ''.join([f"<span style='color:{colors.get(nuc, 'black')}'>{nuc}</span>" for nuc in seq])
    return html

#Read data from a (FASTA) file and create list
def read_file(file_path):
    """
    Read a file and return a list of stripped lines.
    """
    with open(file_path, 'r') as f:
        return ( 
            [line.strip() for line in f] 
        )

def write_file(file_path, seq, mode='w'):
    """
    Write a sequence to a file.
    """
    if isinstance(seq, list):
        seq = '\n'.join(seq)

    with open(file_path, mode) as f:
        f.write(seq + '\n')

    

def parse_fasta(path, seq_type, create_variables=False):
    """
    !! ASSIGN VARIABLE TO FUNCTION TO CREATE LIST !!
    Read a FASTA file with one or more records and returns a list of bio_seq instances, one per record.
    If create_variables is True, create global variables seq_1, seq_2, etc. for each sequence.
    If a sequence is invalid for the given seq_type, it will be skipped with a warning.
    """
    from bio_seq import bio_seq # Avoid circular import
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
        sequences_dict = {}
         #The enumerate function wraps any iterable and yields pairs of (index, element), allowing you to track positions without manual counter management.
        for i, record in enumerate(records): # Loop through each record with its index, index starts at 0 and ends at records length -1
            var_name = f"seq_{i+1}" # Create a variable name like seq_1, seq_2, etc.
            sequences_dict[var_name] = record 
            print(f"Created variable '{var_name}' for sequence with label '{record.label}'")
            print(id(sequences_dict))
        return sequences_dict

    else:
        return records
    
def protein_scrape_uniprot(path=""):
    """Scrape protein sequences from Uniprot given a file with a list of UniProt IDs."""
    with open("", "r") as f:
        temp = [line.strip()[0:6] for line in f]  # Read all lines from the file and strip first 6 characters whitespace/newline characters, you need to call for line in f.readlines() to iterate through each line, else it will treat the whole file as one string
    
    with open("", "r") as f:
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

def parse_fasta_streamlit(file, seq_type, create_variables=False):
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
    import streamlit as st
    from bio_seq import bio_seq
    import io

    lines = file.readlines()
    count = 0
    for line in lines:
        if line.startswith(">"):
            count += 1

    st.write(f"There are {count} total sequences in this file.")

    for line in lines: # Iterate through each line in the file
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
                    st.write(f"Warning: Skipping invalid {seq_type} sequence with label '{label}'")
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
            st.write(f"Warning: Sequence '{label}' is invalid for type '{seq_type}' and will be skipped.")
    st.write(f"Successfully parsed {len(records)} valid {seq_type} sequences from the file.")
    
    if create_variables == True:
        sequences_dict = {}
         #The enumerate function wraps any iterable and yields pairs of (index, element), allowing you to track positions without manual counter management.
        for i, record in enumerate(records): # Loop through each record with its index, index starts at 0 and ends at records length -1
            var_name = f"seq_{i+1}" # Create a variable name like seq_1, seq_2, etc.
            sequences_dict[var_name] = record 
            st.write(f"Created variable '{var_name}' for sequence with label '{record.label}'")
            st.write(id(sequences_dict))
        return sequences_dict

    else:
        return records
    
