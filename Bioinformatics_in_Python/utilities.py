import sys

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
    Read a (FASTA) file and return a list of stripped lines.
    """
    with open(file_path, 'r') as f:
        return ( 
            [line.strip() for line in f] 
        )

def parse_fasta(path, seq_type):
    """
    Read a FASTA file with one or more records.
    Return a list of bio_seq instances, one per record.
    """
    from bio_seq import bio_seq # Avoid circular import
    records = [] #collection of bio_seq objects
    label = None
    seq_parts = [] #list of sequence parts until the next label is found

    with open(path, "r") as fh:
        count = 0
        for line in fh:
            count += line.count('>')
        
        
        for line in fh: # Iterate through each line in the file
            line = line.strip() # Remove leading/trailing whitespace
            if not line: # Skip empty lines
                continue

            if line.startswith(">"): #when a new '>' is found, signals start of a new sequence
                if label is not None: #f label is already set (meaning you just finished reading one record), you:
                    full_seq = "".join(seq_parts) #join all parts of the sequence into one string
                    records.append(bio_seq(full_seq, seq_type, label)) #create a new bio_seq object and add it to the records list
                label = line[1:] #if label is None, set label to the current line without '>'
                seq_parts = [] #reset seq_parts for the next sequence
            else:
                seq_parts.append(line.upper()) # add the current line (part of the sequence) to seq_parts, converting to uppercase

    print(f"There are {count} sequences in this file.")

    if label is not None: #read the last record after the loop ends
        full_seq = "".join(seq_parts)
        records.append(bio_seq(full_seq, seq_type, label))

    return records