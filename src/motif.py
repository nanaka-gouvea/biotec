__author__ = 'natalia'

from origin import generate_neighbourhood
from origin import pool
from origin import pattern_count_aprox
from translation import get_file
from math import log


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
    prob = {n : [0.0 for i in mots[0]] for n in pool}
    for m in mots:
        for i in range(len(m)):
            n = m[i].upper()
            prob[n][i] += 1
    return prob


def profile_motif(mots):
    rows = len(mots)
    profile = count_motif(mots)
    for k, columns in profile.iteritems():
        for i in range(len(columns)):
            profile[k][i] /= rows
    return profile


def consensus(profile):
    columns = profile.values()[0]
    cons = ""
    for i in range(len(columns)):
        max_p = ("", 0.0)
        for k, v in profile.iteritems():
            if v[i] > max_p[1]:
                max_p = (k, v[i])
        cons += max_p[0]
    return cons


def entropy_motif(profile):
    total_entropy = []
    for column in range(len(profile.values()[0])):
        entropy = 0.0
        for k, v in profile.iteritems():
            p = v[column]
            if p > 0:
                entropy += (p * log(p, 2))
        total_entropy.append(abs(entropy))
    return sum(total_entropy)


# motif_mx = [l.replace(" ", "") for l in get_file("/data/motif_matrix.txt").read().splitlines()]
# pmap = profile_motif(motif_mx)
# print entropy_motif(pmap)
