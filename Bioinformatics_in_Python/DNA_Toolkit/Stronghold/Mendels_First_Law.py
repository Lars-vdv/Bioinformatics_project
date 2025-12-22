import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

print(os.getcwd())
file_path = r"C:\Users\larsv\OneDrive\Documenten\VSC\Bioinformatics_project\Bioinformatics_in_Python\DNA_Toolkit\test_data\rosalind_sseq.txt"

from bio_seq import bio_seq
from bio_structs import *
from utilities import *


lst = ['25', '23', '21']
k, m, n = map(float, lst)
t = sum(map(float, lst))
# organize a list with allele one * allele two (possibles) * dominant probability
# multiplications by one were ignored
# remember to substract the haplotype from the total when they're the same for the second haplotype choosed
couples = [
            k*(k-1),  # AA x AA
            k*m,  # AA x Aa
            k*n,  # AA x aa
            m*k,  # Aa x AA
            m*(m-1)*0.75,  # Aa x Aa
            m*n*0.5,  # Aa x aa
            n*k,  # aa x AA
            n*m*0.5,  # aa x Aa
            n*(n-1)*0  # aa x aa
]
# (t-1) indicate that the first haplotype was select
print(round(sum(couples)/t/(t-1), 5))
