__author__ = 'natalia'

from origin import generate_neighbourhood
from origin import pool
from translation import get_file
from origin import pattern_count_aprox


def motif_enumeration(dna, k, d):
    patterns = set()
    nhood = set()
    for seq in dna:
        for i in range(len(seq) - k):
            kmer = seq[i:i + k]
            nhood.update(generate_neighbourhood(pool, k, kmer, d))
    for n in nhood:
        count = 0
        for seq in dna:
            if pattern_count_aprox(seq, n, d):
                count += 1
        if count == len(dna):
            patterns.add(n)
    return patterns