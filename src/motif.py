__author__ = 'natalia'

from origin import generate_neighbourhood
from origin import pool
from origin import pattern_count_aprox
from translation import get_file
from math import log
from origin import hamming_d
from origin import arrange_repeated_dna
from random import randint
from random import seed
from random import randrange
from random import uniform
from random import SystemRandom
from tools.true_random import getnum


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


def profile_motif_pseudo(mots):
    rows = len(mots)
    profile = count_motif(mots)
    for k, columns in profile.iteritems():
        for i in range(len(columns)):
            profile[k][i] = (profile[k][i] + 1) / rows
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


def entropy_column(column):
    e = 0.0
    for p in column:
        if p > 0.0:
            e += (p * log(p, 2))
    return abs(e)


def entropy_motif(profile):
    total_entropy = []
    for column in range(len(profile.values()[0])):
        e = 0.0
        for k, v in profile.iteritems():
            p = v[column]
            if p > 0.0:
                e += (p * log(p, 2))
        total_entropy.append(abs(e))
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


def median_string_many(k, dna):
    min_d = k * len(dna)
    median = []
    possible = arrange_repeated_dna(k)
    for p in possible:
        nd = distance(p, dna)
        if nd < min_d:
            min_d = nd
            median = [p]
        elif nd == min_d:
            median.append(p)
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
    most = ("",0)
    for i in range((len(seq) - k) + 1):
        motif = seq[i:i + k]
        pr = pr_profile(motif, profile)
        if pr > most[1]:
            most = (motif, pr)
    return most[0]


def motifs(profile, dna):
    ms = []
    k = len(profile["A"])
    for seq in dna:
        ms.append(profile_most_probable(seq, k, profile))
    return ms


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
        m_score = score_real_motif(motif, consensus(profile_motif(motif)))
        if m_score < best_motif[1]:
            best_motif = (motif, m_score)
    return best_motif[0]


def greedy_motif_search_pseudo(dna, k):
    t = len(dna)
    best_motif = ([], t * k)
    for seq in dna:
        best_motif[0].append(seq[0:k])
    for i in range(len(dna[0]) - k):
        motif = [dna[0][i:i + k]]
        for seq in dna[1:]:
            profile = profile_motif_pseudo(motif)
            kmer = profile_most_probable(seq, k, profile)
            motif.append(kmer)
        m_score = score_real_motif(motif, consensus(profile_motif_pseudo(motif)))
        if m_score < best_motif[1]:
            best_motif = (motif, m_score)
    return best_motif[0]


def random_bias(bias_list):
    number = uniform(0, sum(bias_list))
    current = 0
    for i, bias in enumerate(bias_list):
        current += bias
        if number <= current:
            return i


def randomized_motif_search(dna, k):
    random_motifs = []
    for seq in dna:
        si = randint(1,len(dna[0]) - k)
        random_motifs.append(seq[si:si + k])
    best_motif = (random_motifs, score_real_motif(random_motifs, consensus(profile_motif_pseudo(random_motifs))))
    while True:
        prof = profile_motif_pseudo(random_motifs)
        random_motifs = motifs(prof, dna)
        cons = consensus(prof)
        sc = score_real_motif(random_motifs, cons)
        if sc < best_motif[1]:
            best_motif = (random_motifs, sc)
        else:
            return best_motif


def randomized_motif_search_times(k, t, dna, n):
    seed(8)
    best = ([], t * k)
    for i in range(n):
        n_best = randomized_motif_search(dna, k)
        if n_best[1] <= best[1]:
            best = (n_best[0], n_best[1])
    return best


def gibs_sampling(k, t, n, dna):
    random_motifs = []
    seq_len = len(dna[0])
    for seq in dna:
        si = randrange(seq_len - k + 1)
        random_motifs.append(seq[si:si + k])
    best_motif = (random_motifs, score_real_motif(random_motifs, consensus(profile_motif_pseudo(random_motifs))))
    for _ in xrange(n):
        i = randrange(t)
        random_motifs.pop(i)
        prof = profile_motif_pseudo(random_motifs)
        ith_seq = dna[i]
        probs = [pr_profile(ith_seq[s:s + k], prof) for s in range(seq_len - k + 1)]
        ex = random_bias(probs)
        motif_i = ith_seq[ex:ex + k]
        random_motifs.insert(i, motif_i)
        #TODO usar dna em vez de motif para encontrar o consensus?
        sc = score_real_motif(random_motifs, consensus(profile_motif_pseudo(random_motifs)))
        if sc < best_motif[1]:
            best_motif = (random_motifs, sc)
    return best_motif


def gibs_sampling_times(k, t, n, m, dna):
    best = ([], t * k)
    for i in xrange(m):
        # seed(i)
        n_best = gibs_sampling(k, t, n, dna)
        if n_best[1] < best[1]:
            best = (n_best[0], n_best[1])
    return best

# seed(31)
# dnas = "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA".split(" ")
# print ' \n'.join(gibs_sampling_times(8, 5, 100, 20, dnas)[0])

# expected = "TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG".split(" ")
# for r in range(1000):
#     seed(r)
#     if set(gibs_sampling_times(8, 5, 100, 20, dnas)[0]) == set(expected):
#         print r
#         break
# dnas = get_file("/data/gibs_test.txt").read().splitlines()
# print ' \n'.join(gibs_sampling_times(15, 20, 5000, 20, dnas)[0])

