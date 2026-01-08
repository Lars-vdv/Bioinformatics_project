import itertools

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

print(os.getcwd())
file_path = r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\test_data\rosalind_tran.txt"

from bio_seq import bio_seq
from bio_structs import *
from utilities import *

def enumerater(number):
    lst = []
    for i in range(1, number + 1):
        lst.append(i)
    permutations = list(itertools.permutations(lst))
    print(len(permutations))
    for i in range(len(permutations)):
        permutations[i] = list(permutations[i])
        for j in range(len(permutations[i])):
            permutations[i][j] = str(permutations[i][j])
    for i in permutations:
        print(' '.join(i))

test = enumerater(7)


write_file("test.txt", test, mode='w')