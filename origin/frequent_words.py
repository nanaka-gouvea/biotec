# coding=utf-8
from tools import OrderedSet
import collections

__author__ = 'natalia'


complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}


def pattern_count(text, pattern):
    count = 0
    for i in range(len(text) - (len(pattern) - 1)):
        if text[i:i + len(pattern)] == pattern:
            count += 1
    return count


def most_frequent_words(text, size):
    frequency = {}
    max_count = 0
    frequents = []
    for i in range(len(text) - (size - 1)):
        word = text[i:i + size]
        frequency[word] = frequency[word] + 1 if word in frequency else 1
        count = frequency[word]
        if count > max_count:
            max_count = count
            frequents = [word]
        else:
            if count == max_count:
                frequents.append(word)
    return ' '.join(frequents)


def reverse_complement(strand):
    rc = []
    i = len(strand)
    while i > 0:
        i -= 1
        rc.append(complement[strand[i]])
    return ''.join(rc)


def pattern_starts(pattern, text):
    starts = []
    for i in range(len(text)-(len(pattern)-1)):
        if text[i:i + len(pattern)] == pattern:
            starts.append(i)
    return ' '.join(str(x) for x in starts)

def skew_g_c(genome):
    skew = [0]
    i = 0
    while i < len(genome):
        update = skew[i]
        if genome[i] == 'G':
            update += 1
        if genome[i] == 'C':
            update -= 1
        skew.append(update)
        i += 1
    return skew


def minimum_skew(genome):
    prev = 0
    minimum = 0
    min_loc = []
    i = 0
    while i < len(genome):
        # update = prev
        current = genome[i]
        if current == 'G':
            prev += 1
        if current == 'C':
            prev -= 1
        # skew.append(update)
        if prev < minimum:
            minimum = prev
            min_loc = [i + 1]
        elif prev == minimum:
            min_loc.append(i + 1)
        i += 1
    return min_loc

def maximum_skew(genome):
    prev = 0
    maximum = 0
    max_loc = []
    i = 0
    while i < len(genome):
        current = genome[i]
        if current == 'G':
            prev += 1
        if current == 'C':
            prev -= 1
        if prev > maximum:
            maximum = prev
            max_loc = [i + 1]
        elif prev == maximum:
            max_loc.append(i + 1)
        i += 1
    return max_loc


def hamming(a, b):
    i = 0
    count = 0
    while i < len(a):
        if a[i] != b[i]:
            count += 1
        i += 1
    return count


def pattern_starts_aprox(pattern, text, d):
    starts = []
    for i in range(len(text)-(len(pattern)-1)):
        word = text[i:i + len(pattern)]
        if hamming(word, pattern) <= d:
            starts.append(i)
    return ' '.join(str(x) for x in starts)


def pattern_count_aprox(text, pattern, d):
    count = 0
    for i in range(len(text) - (len(pattern) - 1)):
        word = text[i:i + len(pattern)]
        if hamming(word, pattern) <= d:
            count += 1
    return count


def most_frequent_words_aprox(text, size, d):
    frequency = {}
    max_count = 0
    frequents = []
    for i in range(len(text) - (size - 1)):
        word = text[i:i + size]
        if word not in frequency:
            frequency[word] = 0
        matches = False
        for w in frequency.keys():
            if hamming(w, word) <= d:
                matches = True
        if matches:
            frequency[word] += 1
            count = frequency[word]
            if count > max_count:
                max_count = count
                frequents = [word]
            else:
                if count == max_count:
                    frequents.append(word)
    return ' '.join(frequents)

def anagrams(word):
    if len(word) <= 1:
        return word
    result = []
    for ana in anagrams(word[1:]):
        for i in range(len(word)):
            result.append(ana[:i] + word[0] + ana[i:])
    return sorted(result)

def arrange_repeated(letters, size):
    if size == 0:
        return []
    result = []
    if size == 1:
        for l in letters:
            result.append(l)
        return result
    for a in arrange_repeated(letters, size - 1):
        for i in range(len(letters)):
            result.append(letters[i] + a[:i] + a[i:])
    return sorted(result)

pool = "ACGT"
def pattern_to_number(pattern):
    possible = arrange_repeated(pool, len(pattern))
    for i in range(len(possible)):
        if possible[i] == pattern:
            return i
# print pattern_to_number("ATGCAA")

def number_to_pattern(i, size):
    possible = arrange_repeated(pool, size)
    return possible[i]

# print number_to_pattern(5437, 8)

def frequency_array(seq, size):
    #overlapping counts!
    possible = arrange_repeated(pool, size)
    fa = collections.OrderedDict()
    for p in possible:
        fa[p] = 0
    for i in range(len(seq) - (size - 1)):
        fa[seq[i:i + size]] += 1
    return fa.values()


def clump_finding(genome, k, l, t):
    # fa = {word:[starting indexes]}
    fa = {}
    clump_patterns = set()
    i = k
    while i <= len(genome):
        si = i - k
        word = genome[si:i]
        try:
            fa[word].append(si)
        except KeyError:
            fa[word] = [si]
        fq = fa[word]
        if len(fq) >= t:
            #start index of the first match in the last t matches
            th_match = fq[len(fq) - t]
            if i - th_match <= l:
                clump_patterns.add(word)
        i += 1
    return clump_patterns

def clump_finding_not_overlapping(genome, k, l, t):
    # fa = {word:[starting indexes]}
    fa = {}
    clump_patterns = set()
    i = k
    while i <= len(genome):
        si = i - k
        word = genome[si:i]
        try:
            fa[word].append(si)
        except KeyError:
            fa[word] = [si]
        fq = fa[word]
        len_fq = len(fq)
        if len_fq >= t:
            j = len_fq - 1
            lsi = fq[j]
            matches = 1
            while j > 0:
                if fq[j - 1] < i - l:
                    break
                if lsi - fq[j - 1] >= k:
                    matches += 1
                    lsi = fq[j - 1]
                j -= 1
            if matches >= t:
                clump_patterns.add(word)
                #start index of the first match in the last t matches
                # th_match = fq[len_fq - t]
                # if i - th_match <= l:
                #     overlapping = False
                #     for ix in range(t + 1)[1:]:
                #         if fq[-ix] - fq[-ix - 1] > k:
                #             overlapping = True
                #     if not overlapping:
                #         clump_patterns.add(word)
        i += 1
    return clump_patterns


def smart_neighbourhood_count(seq, max_hamming):
    size = len(seq)
    return pow(size, size) - pow(max_hamming, size)

def neighbourhood_count(seq, max_hamming):
    count = 0
    for i in arrange_repeated(pool, len(seq)):
        if hamming(i, seq) <= max_hamming:
            count += 1
    return count

print neighbourhood_count("ATTA", 3)
print smart_neighbourhood_count("ATTA", 3)
# print len(clump_finding_not_overlapping(open("data/E-coli.txt").read(), 9, 500, 3))
#1472
# clump_finding_not_overlapping(open("data/E-coli.txt").read(), 9, 500, 3)
# print len(find_clump_patterns(open("data/E-coli.txt").read(), 10, 200, 4))
# print ' '.join(new_find_clump_patterns(open("data/E-coli.txt").read(), 10, 200, 4))
# print len(find_clump_patterns(open("data/E-coli.txt").read(), 9, 500, 4))
# print len(find_clump_patterns(open("data/E-coli.txt").read(), 9, 500, 5))
# print len(find_clump_patterns(open("data/E-coli.txt").read(), 9, 500, 6))
# For t=4 I get 588;for t=5 I get 288;t=6 I get 34;

# other = "AATGATGAAA TGGGTCAAAA GGGTCAAAAG ATGATGAAAT AGCCGCAACA CGCAACAACC CCGCAACAAC GGTCAAAAGT GAGTTAAATA TTAAATAATC TGAAATGATG GCAGCCGCAA ACTATGGCAC GTTAAATAAT AAAGTTGCCG TTATCCCCGC CCGCTGGCGC AAGAAAGCGG GACGGTGCTA CGGGGAACTC ATCCCCGCTG GCGCGGGGAA CACTATGGCA ATCAGCAGCC CCCGCTGGCG GGCGCGGGGA TATGGCACTA TGGCGCGGGG GCGGGGAACT GTCAAAAGTT AGTTAAATAA CGGTTTATCC GAAATGATGA GGCACTATGG TCAAAAGTTG GTGGGTCAAA CAGCCGCAAC TATCCCCGCT CGGGGAACAC CTATGGCACT TATCAGCAGC TGGCACTATG GGTGGGTCAA GATGAAATGA TGATGAAATG TTTATCCCCG AGGAGTTAAA GGAGTTAAAT CAGCAGCCGC GCTGGCGCGG AAATGATGAA GCCGCAACAA ATGAAATGAT AGCAGCCGCA AAAAGTTGCC CGCTGGCGCG GGACGGTGCT GCACTATGGC CTGGCGCGGG GTTTATCCCC CAGGAGTTAA GCGGGGAACA ATGGCACTAT CCCCGCTGGC CGCGGGGAAC CAAAAGTTGC TCCCCGCTGG GGTTTATCCC".split(' ')
# mine =  "CGCAACAACC GGTTTATCCC GACGGTGCTA ATGATGAAAT TTTATCCCCG GCGGGGAACT GATGAAATGA TCAAAAGTTG CGGGGAACTC AGCAGCCGCA GGTGGGTCAA AAATGATGAA TGATGAAATG CAGCCGCAAC ATCCCCGCTG CGCGGGGAAC GGCGCGGGGA GCTGGCGCGG CAGCAGCCGC AAAAGTTGCC TTATCCCCGC CGGTTTATCC CCCGCTGGCG GCACTATGGC ATCAGCAGCC CACTATGGCA GCGGGGAACA CTGGCGCGGG AGCCGCAACA GTTTATCCCC AATGATGAAA ATGAAATGAT CCGCAACAAC TGGCACTATG TATGGCACTA CCCCGCTGGC CTATGGCACT GCAGCCGCAA ACTATGGCAC CAAAAGTTGC GGGTCAAAAG CGCTGGCGCG GGCACTATGG GAAATGATGA CCGCTGGCGC GCCGCAACAA CGGGGAACAC GCGCGGGGAA TGGCGCGGGG TCCCCGCTGG TGGGTCAAAA GTGGGTCAAA GGACGGTGCT GTCAAAAGTT GGTCAAAAGT AAAGTTGCCG TATCAGCAGC ATGGCACTAT TATCCCCGCT TGAAATGATG AAGAAAGCGG".split(' ')
# print ' '.join(set(mine).difference(other))
# print ' '.join(set(mine).intersection(other))
# common = set(set(mine).intersection(other))
#
# for x in other:
#     if x not in mine:
#         print x

#didnt find these clumps:
# GAGTTAAATA
# TTAAATAATC
# GTTAAATAAT
# AGTTAAATAA
# AGGAGTTAAA
# GGAGTTAAAT
# CAGGAGTTAA


#OK! all clumps i found really are clumps!
# for x in mine:
#     if x not in other:
#         print x
