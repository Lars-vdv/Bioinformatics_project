import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

print(os.getcwd())
file_path = r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\test_data\rosalind_gbk.txt"


from bio_seq import bio_seq
from bio_structs import *
from utilities import *

file = read_file(file_path) 

genus_name = str(file[0])
start_date = str(file[1])
end_date = 0

import webbrowser

#use triple quotes to avoid issues with special characters in the URL
# test_url = """www.ncbi.nlm.nih.gov/nuccore?term="Entosthodon"[Organism] AND ("2004/08/02"[PDAT] : "2012/12/03"[PDAT])&cmd=DetailsSearch"""
# webbrowser.open(test_url, new=2, autoraise=True)

def search_genbank(organism, start_date = 0, end_date = 0):
    """Searches and opens GenBank for the given organism and optional date range."""

    if start_date != 0 and end_date != 0:
        url = f"""www.ncbi.nlm.nih.gov/nuccore?term="{organism}"[Organism] AND ("{start_date}"[PDAT] : "{end_date}"[PDAT])&cmd=DetailsSearch"""

    if start_date == 0 and end_date == 0:
        url = f"""www.ncbi.nlm.nih.gov/nuccore?term="{organism}"[Organism]&cmd=DetailsSearch"""

    if start_date != 0 and end_date == 0:
        url = f"""www.ncbi.nlm.nih.gov/nuccore?term="{organism}"[Organism] AND ("{start_date}"[PDAT] : "3000"[PDAT])&cmd=DetailsSearch"""

    if start_date == 0 and end_date != 0:
        url = f"""www.ncbi.nlm.nih.gov/nuccore?term="{organism}"[Organism] AND ("0000"[PDAT] : "{end_date}"[PDAT])&cmd=DetailsSearch"""

    webbrowser.open(url, new=2, autoraise=True)

#search_genbank(genus_name, start_date, end_date)

# Function to scrape GenBank for accession numbers
def scrape_genbank_accession(accession_number="", file_path="", email="your_name@your_mail_server.com"):
    """Scrapes GenBank for the given accession number and returns FASTA.txt file."""
    """Enter accession_number or file_path to a text file with accession numbers (one per line)."""

    from Bio import SeqIO, Entrez

    if accession_number == "" and file_path == "":
        print("Please provide either an accession number or a file path.")
        return
    
    if accession_number != "" and file_path != "":
        print("Please provide either an accession number or a file path, not both.")
        return
    
    if accession_number == "" and file_path != "":
        with open(file_path, "r") as f:
            temp = [line.strip() for line in f]

    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=temp, rettype="fasta")
    records = handle.read()
    write_file("GenBank_FASTA.txt", records)
    
scrape_genbank_accession(accession_number="", file_path=r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\Stronghold\test.txt", email="larsvdvelden@live.nl")

# Parsing the downloaded FASTA file and printing sequence info
seqs = bio_seq()
seqs = parse_fasta("GenBank_FASTA.txt", seq_type="DNA", create_variables=True)

# for i in range(len(seqs.keys())):
    
#     print(bio_seq.get_seq_info(seqs.get(f"seq_{i+1}")))

# find shortest sequence for Rosalind problem

lst = []
for i in range(len(seqs.keys())):
    x = bio_seq.get_seq_len(seqs.get(f"seq_{i+1}"))
    lst.append(x)

min_val = min(lst)
min_idx, min_val = min(enumerate(lst), key=lambda iv: iv[1])
#enumerate takes sequence and returns index-value pairs
# key parameter specifies to use value for comparison,tells min() how to compare items
# lambda iv: iv[1] is an anonymous function that takes an index-value pair and returns the 
# value for comparison
# lambda functions are used for short, throwaway functions that are not reused elsewhere 
# so in this case it is similar to defining a function like:
# def get_value(iv):
#     return iv[1]

print(min_val, min_idx)

output = bio_seq.get_seq_info(seqs.get(f"seq_{min_idx+1}"))
print(output)

# obj = seqs.get(f"seq_{min_idx+1}")
# label = obj.label
# sequence = obj.seq

test = bio_seq(seq="ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", label="Test_Sequence_1")

# write_file("Stronghold_GenBank_output.txt", f">{label}\n{sequence}\n")

write_fasta("Stronghold_GenBank_output.fasta", label = test.label,sequence= test.seq)