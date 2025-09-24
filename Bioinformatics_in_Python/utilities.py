import sys

def colored(seq):
    """ Change ATCG to a colored letter """
    colors = {
    'A': '\033[31m',  # red
    'C': '\033[12m',  # blue
    'G': '\033[33m',  # yellow
    'T': '\033[32m',  # green
    'U': '\033[35m',  # magenta
    'reset': '\033[0m'
    }

    result = ""
    for nuc in seq:
        color_code = colors.get(nuc, colors['reset'])
        result += f"{color_code}{nuc}{colors['reset']}"
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
