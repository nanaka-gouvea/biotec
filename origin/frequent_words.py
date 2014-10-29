__author__ = 'natalia'
import os
from os import path
import collections
# coding=utf-8

complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}


def pattern_count(text, pattern):
    count = 0
    for i in range(len(text)-(len(pattern)-1)):
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
            min_loc = [i+1]
        elif prev == minimum:
            min_loc.append(i+1)
        i += 1
    return min_loc


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

#working for all except large and forum
def old_working_find_clump_patterns(genome, k, l, t):
    #frequency = word:((laststartindex, lastendindex), count in this window)
    frequency = {}
    clump_patterns = set()
    i = k
    while i <= len(genome): # vai ate o 'len', pois o i sera usado para intervalo exclusivo
        si = i-k #start index
        word = genome[si:i]
        count = 1
        try:
            f = frequency[word]
        except KeyError:
            frequency[word] = ((si, i), count)
        else:
            # if the word's last start index is inside my current window, increment counter
            #(endindex is exclusive)
            if f[0][0] >= i - l: #in window
                count = f[1] + 1
                frequency[word] = ((f[0][0], i), count)
            else: #not in window, restart counting
                frequency[word] = ((si, i), count)
        if count == t: #if counter in this window reaches minimun, it's a winner
            clump_patterns.add(word)
        i += 1
    return clump_patterns

def old_not_working_find_clump_patterns(genome, k, l, t):
    # if my LAST si is not in window, I MUST start counting again (and reset indexes)
    # if my first starting index is not in window, pop first si, dont increment count

    #frequency = word:(starting indexes, count in this window, starting indexes index, last si)
    frequency = {}
    clump_patterns = set()
    i = k
    # vai ate o 'len', pois o i sera usado para intervalo exclusivo
    while i <= len(genome):
        #start index
        si = i - k
        word = genome[si:i]
        count = 1
        try:
            f = frequency[word]
        except KeyError:
            frequency[word] = [[si], count, 0, si]
        else:
            # update start indexes
            start_indexes = f[0]
            start_indexes.append(si)

            count = f[1]
            #f[2] = start index index, uses index as parameter instead of 'pop'ing firsts, for performance
            f_start = start_indexes[f[2]]
            #if last start index not in window, restart counting
            if f[3] <= i - l:
                count = 1
                frequency[word] = [[si], count, 0, si]
            else:
                # if the word's first start index is inside my current window, increment count ((endindex is exclusive))
                if f_start >= i - l:
                    #in window
                    count += 1
                    frequency[word][1] = count
                else:
                    #not in window, do not alter counting, update first start index (first start index moves to next)
                    frequency[word][2] += 1
        if count == t:
            #if counter in this window reaches minimun, it's a winner
            clump_patterns.add(word)
        i += 1
    return clump_patterns

def find_clump_patterns(genome, k, l, t):
    # pkmers = arrange_repeated(pool, k)
    # fa = {p:[] for p in pkmers}

    fa = {}
    clump_patterns = set()
    i = k
    # for i in range(len(genome) - (k - 1)):
    while i <= len(genome):
        si = i - k
        word = genome[si:i]
        try:
            fa[word].append(si)
        except KeyError:
            fa[word] = [si]
        else:
            fq = fa[word]
            if len(fq) >= t:
                last_fq_ix = len(fq) - 1
                #TODO melhorar essa verificacao, quando ainda nao atingiu a janela t se torna negativo
                todo = fq[last_fq_ix - t]
                primeiro_t_ix = todo if todo >= 0 else 0
                if i - primeiro_t_ix <= l:
                    clump_patterns.add(word)
        i += 1
    return clump_patterns

# print len(arrange_repeated("ACGT", 9))

# print find_clump_patterns("GCTGGCTGGCTGGCTGG", 9, 17, 3)
# print find_clump_patterns("AGAGAGAGAGA", 3, 9, 4)
# print find_clump_patterns("AGAGA", 3, 4, 2)
# print find_clump_patterns("AGAG", 1, 4, 2)
print find_clump_patterns("AGAG", 2, 3, 2)
# print find_clump_patterns("TAG", 3, 3, 1)
# print os.chdir("data")
# print len(find_clump_patterns(open("data/test_clump.txt").read(), 12, 595, 19))
# print ' '.join(find_clump_patterns(open("data/test_clump.txt").read(), 12, 595, 19))
# AGAGTGATTGCG GTGGATAGCCTA GTGATCCACCGA GATAGTTGGTCT ACTTCCAAACAG TACTCCTGAAGT TTGCAAACTGAC CCGCACGAAGTA ATAACGATTTCC
# print ' '.join(find_clump_patterns("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", 5, 50, 4))

# print len(find_clump_patterns(open("data/E-coli.txt").read(), 9, 500, 3))
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
