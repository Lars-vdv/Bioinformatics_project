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
test_url = """www.ncbi.nlm.nih.gov/nuccore?term="Entosthodon"[Organism] AND ("2004/08/02"[PDAT] : "2012/12/03"[PDAT])&cmd=DetailsSearch"""
webbrowser.open(test_url, new=2, autoraise=True)

def parse_genbank(organism, start_date = 0, end_date = 0):
    
    if start_date != 0 and end_date != 0:
        url = f"""www.ncbi.nlm.nih.gov/nuccore?term="{organism}"[Organism] AND ("{start_date}"[PDAT] : "{end_date}"[PDAT])&cmd=DetailsSearch"""

    if start_date == 0 and end_date == 0:
        url = f"""www.ncbi.nlm.nih.gov/nuccore?term="{organism}"[Organism]&cmd=DetailsSearch"""

    if start_date != 0 and end_date == 0:
        url = f"""www.ncbi.nlm.nih.gov/nuccore?term="{organism}"[Organism] AND ("{start_date}"[PDAT] : "3000"[PDAT])&cmd=DetailsSearch"""

    if start_date == 0 and end_date != 0:
        url = f"""www.ncbi.nlm.nih.gov/nuccore?term="{organism}"[Organism] AND ("0000"[PDAT] : "{end_date}"[PDAT])&cmd=DetailsSearch"""

    webbrowser.open(url, new=2, autoraise=True)

parse_genbank(genus_name, start_date, end_date)
