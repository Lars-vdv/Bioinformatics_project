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
