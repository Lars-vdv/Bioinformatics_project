#!/usr/bin/env python3
import argparse
from typing import NamedTuple, TextIO
from Bio import SeqIO


class Args(NamedTuple):
    file: TextIO


def get_args() -> Args:
    parser = argparse.ArgumentParser(description='longest substring', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('rosalinf_lcsm.txt', metavar='FASTA FILE', type=argparse.FileType('rt'), help=r'test_data/rosalind_lcsm.txt"')

    args = parser.parse_args()

    return Args(args.file)


def main():
    args = get_args()
    seqs = SeqIO.parse(args.file, 'fasta')
    seqs_list = [str(rec.seq) for rec in seqs]
    min_len = len(min(seqs_list))

    shared_kmers = set()
    for k in range(min_len, 0, -1):
        list_of_kmers = []
        with open(args.file.name, 'rt') as file:
            seqs = SeqIO.parse(file, 'fasta')
            for rec in seqs:
                list_of_kmers.append(find_kmers(str(rec.seq), k))
        shared = set.intersection(*map(set, list_of_kmers))
        if shared:
            shared_kmers.update(shared)
            break

    print(shared_kmers)


def find_kmers(seqs, k: int):
    n = len(seqs) - k + 1
    return [] if n < 1 else [seqs[i:i + k] for i in range(n)]


if __name__ == '__main__':
    main()
