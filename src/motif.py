__author__ = 'natalia'

from origin import generate_neighbourhood
from origin import pool
from origin import pattern_count_aprox
from translation import get_file
from math import log
from origin import hamming_d
from origin import arrange_repeated_dna


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
    prob = {n: [0.0 for i in mots[0]] for n in pool}
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
            if p > 0.0:
                entropy += (p * log(p, 2))
        total_entropy.append(abs(entropy))
    return sum(total_entropy)


def distance(pattern, dna):
    k = len(pattern)
    d = 0
    for s in dna:
        hd = k
        for i in range(len(s) - k):
            p = s[i:i + k]
            nhd = hamming_d(pattern, p)
            if nhd < hd:
                hd = nhd
        d += hd
    return d


def median_string(k, dna):
    min_d = k * len(dna)
    median = ""
    possible = arrange_repeated_dna(k)
    for p in possible:
        nd = distance(p, dna)
        if nd < min_d:
            min_d = nd
            median = p
    return median


def score_real_motif(mots, cons):
    count = 0.0
    for m in mots:
        count += hamming_d(m, cons)
    return count


def pr_profile(seq, profile):
    pr = 1.0
    for i in range(len(seq)):
        pr *= profile[seq[i]][i]
        if pr == 0.0:
            return pr
    return pr


def profile_most_probable(seq, k, profile):
    most = (seq[0:k],0)
    for i in range(len(seq) - k):
        motif = seq[i:i + k]
        pr = pr_profile(motif, profile)
        if pr > most[1]:
            most = (motif, pr)
    return most[0]


def greedy_motif_search(dna, k):
    t = len(dna)
    best_motif = ([], t * k)
    for seq in dna:
        best_motif[0].append(seq[0:k])
    for i in range(len(dna[0]) - k):
        motif = [dna[0][i:i + k]]
        for seq in dna[1:]:
            profile = profile_motif(motif)
            kmer = profile_most_probable(seq, k, profile)
            motif.append(kmer)
        # m_score = score_motif(motif)
        m_score = score_real_motif(motif, consensus(profile_motif(motif)))
        if m_score < best_motif[1]:
            best_motif = (motif, m_score)
    return best_motif[0]

# motif_mx = [l.replace(" ", "") for l in get_file("/data/motif_matrix.txt").read().splitlines()]
# pmap = profile_motif(motif_mx)
#
# for k, v in pmap.iteritems():
#     print k, ": ", v
#
# print pr_profile("TCGTGGATTTCC", pmap)

# pcolumns = [[float(x) for x in l.split(" ")] for l in get_file("/data/p_matrix_test.txt").read().splitlines()]
# test_p = {"A":pcolumns[0], "C":pcolumns[1], "G":pcolumns[2], "T":pcolumns[3]}
# print profile_most_probable("TCGATCCGAGCGCCATATGCTATAAGTCGCTAACGTTTAGATAGCCCGCTGCGAATCAACGTTAAGATATTCATACCTGGTATACTCCGGAGTAGGTTGTGTTAGCACCCTTTTGAGCGACGTCACGTGCGAGTAGGGCTGGGGACCGTAAACCAGGCTGTCTCATTGACACAGCCAATAGTAATTAGGCCTTCGAGACTATCTCCGTTTTCATGTGCTCCCAACGCATTAGTCTCCATCAGACGATGCTTCCGGTTGTAGGCCCTGAACCCAGGTTCGCCCCATACGACAATGCGTTCTCCCGTCCTTTTTCCCAAGTTTGATCTCAGCCGCTCGGCTTACGCGTAACAAGTACACGACATTCCATCACCGCCAAGTGAAGAACGCATTTATTCTCTCCACGTAGAACCAATATTCTGCGATCTAACGACAAAACGCACAAGTTCGACGCCAGTCGGCCATACTAGTCGAATCGTCCTCTTCAAAGTGACATTTTCGACCGCAAATTTCGTGCCCCTGCAGCCTTACAAGATCCGTATTCCGATCTTATCTGTGGCGGCCATGGGGTTCCCACGTCTCCGCGTAACATGTCCAAATTCAGCCTACATTTGGCTAGGGTATGGGATCATGGGTTGTCTAAATTGTGGATAGCAAACCAAGTGTATGAAACGTTATAGACACGTGCCGCATACACCGAGATGACATGCAGCTCGAACGGCCATGCTTACGGGGCGGTGGAAGCTTCGAATCCGCACGCCGACCCGTGCTGCGAATTGTACCAAAGAGGTGTGACCCGCACCTAATGGGTTACGGACTCCAGCTACACTCTCATGTGGGTGCCAGGCGGGCAGTTCTCGGCGGCTCCAGATTACAGCAGGAATCGGACGCATTGAAGACTTTACCGCGTTGCCCATCAACGAGTGCCCCTGATGTTACCAGCTAAAAGCACAGGAATCTTTGCCGTCCGTGGTCGTAGTAATACGCCATCCGGCCACCAAAG", 15, test_p)
# print median_string(6, "ACGAGGGATAGTAATGCCTTGACGATGATATAGGAAGCGCGT AGACATATCATCCTCAGGTCCTAACAATAGGCGAGGTATGCA TCGAGGTACGGGTACACAGCTTGGATGCGGTGGGATGAAGGA CCGAGGGTCATTGCTGCGAACGTTTGACGCGTCTCTGTTCCT TGGATGAAGAATTATGAGATTCGACGTCGGGCGAGGAGTTAT TAGCTGAGGAACTTGACGTCCAGCTCGAGGCGTGGGCCCTTG GGGATCACGAGGGCGGCAAACAGCAAAAATCTTCCGATTACA TGACTCACATGTTCTCGCCTTTTATCAAAGCCGAGGTTCCTT GCGAGGACCAAAGAGCTGGCCGTTCGGTCTGGTTCGTCTATA CGAAGAATAGTGAAAATATTTTGTGCGAGGGTCCGGTTAGAG".split(" "))
# GCGAGG

dnas = get_file("/data/dna_test.txt").read().splitlines()
print ' \n'.join(greedy_motif_search(dnas, 12))
# result = ['AGTGGGTATCTC', 'TAAAAAGGTATA', 'AACCACGAGTAC', 'TGTCATGTGCGG', 'AACCTAAACCCT', 'AGTCGTTATCCC', 'AGTAATATGTAC', 'AGTGGTTATCAC', 'AGTGGTTATCCC', 'AGTGGCTATCGC', 'AGTGGATATCCC', 'AGTGAGAAGCAA', 'AGTGACTAGACA', 'TAAGACTAGTTA', 'TATGAAGGGTGA', 'AGTCGGGATAAC', 'AGTGGGTATCTC', 'AGCGGTTAGTCA', 'AGTGAAATTCCT', 'TGTGGATGGCTT', 'TGTAGGTATCAC', 'TGCAGATATCCA', 'TGTGGTTATCAC', 'TGTCATTATTCA', 'TGCGTAGATCAA']
# expected = "AGTGGGTATCTC TAAAAAGGTATA AACCACGAGTAC TGTCATGTGCGG AACCTAAACCCT AGTCGTTATCCC AGTAATATGTAC AGTGGTTATCAC AGTGGTTATCCC AGTGGCTATCGC AGTGGATATCCC AGTGAGAAGCAA AGTGACTAGACA TAAGACTAGTTA TATGAAGGGTGA AGTCGGGATAAC AGTGGGTATCTC AGCGGTTAGTCA AGTGAAATTCCT TGTGGATGGCTT TGTAGGTATCAC TGCAGATATCCA TGTGGTTATCAC TGTCATTATTCA TGCGTAGATCAA".split(" ")
