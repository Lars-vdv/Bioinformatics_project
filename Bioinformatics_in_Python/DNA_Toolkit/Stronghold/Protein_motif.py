import os
print(os.getcwd())  # Check current working directory

from bio_seq import *
from bio_structs import *
from utilities import *
import urllib.request

seq_1 = urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/B5ZC00.fasta")
print(type(seq_1))