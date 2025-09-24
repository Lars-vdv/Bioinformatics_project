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


#Read data from a FASTA file and create list
def read_file(file_path):
    """Reading a (FASTA) file and returning a list of lines"""
    with open(file_path, 'r') as f: 
        # 'r' opens the file for reading. 'r', 'w' writing, truncating the file first 'x' create new file and open for writing 'a' open for writing, appending to the end of the file if it exists
        # 'as f' assigns the open file object to f
        return [l.strip() for l in f.readlines()]
        # f.readlines() returns a list where each element is one line of the file, including the trailing newline \n
        #Iterates over each raw line l from f.readlines(), Calls l.strip(), which removes leading and trailing whitespace (including \n)


#Read data from a FASTA file and create list
def read_file(file_path):
    """Reading a (FASTA) file and returning a list of lines"""
    with open(file_path, 'r') as f: 
        # 'r' opens the file for reading. 'r', 'w' writing, truncating the file first 'x' create new file and open for writing 'a' open for writing, appending to the end of the file if it exists
        # 'as f' assigns the open file object to f
        return [l.strip() for l in f.readlines()]
        # f.readlines() returns a list where each element is one line of the file, including the trailing newline \n
        #Iterates over each raw line l from f.readlines(), Calls l.strip(), which removes leading and trailing whitespace (including \n)
