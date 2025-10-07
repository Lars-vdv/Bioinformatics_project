nucleotides = ["A", "C", "G", "T"] #[] make list
reverse_dna_seq = {'A' : 'T', 'T' : "A", 'G' : 'C', 'C' : 'G'} 
DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

MW_proteins = {"A": 71.037113805, "C": 103.009184505, "D": 115.026943505,
               "E": 129.042593505, "F": 147.068413945, "G": 57.021463505,
               "H": 137.058911875, "I": 113.084064015, "K": 128.094963050,
                "L": 113.084064015, "M": 131.040484645, "N": 114.042927470,
                "P": 97.052763875, "Q": 128.058577540, "R": 156.101111050,
                "S": 87.032028435, "T": 101.047678505, "V": 99.068413945,
                "W": 186.079312980, "Y": 163.063328575, "_": 0.0}
