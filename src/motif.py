__author__ = 'natalia'

from origin import generate_neighbourhood
from origin import pool
from origin import pattern_count_aprox
from translation import get_file


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


lower_pool = "acgt"


def score_motif(mots):
    count = 0
    for m in mots:
        for l in m:
            if l in lower_pool:
                count += 1
    return count

def count_motif(mots):
    prob = {"A": {}, "C":{}, "G":{}, "T":{}}
    for m in mots:
        for i in range(len(m)):
            n = m[i].upper()
            try:
                prob[n][i] += 1
            except KeyError:
                prob[n][i] = 1
    return prob

motif_mx = [l.replace(" ", "") for l in get_file("/data/motif_matrix.txt").read().splitlines()]
cmap = count_motif(motif_mx)

for k, v in cmap.iteritems():
    print k, ": ", ' '.join([str(x) for x in v.values()])


