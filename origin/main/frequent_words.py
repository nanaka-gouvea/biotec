__author__ = 'natalia'
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

def find_clump_patterns(genome, k, l, t):
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
            # if the word's last start index is inside my current window, and the last end index is before my current start (no overlapping), increment counter
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

def new_find_clump_patterns(genome, k, l, t):
    #frequency = word:(starting indexes, count in this window, starting indexes index)
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
            frequency[word] = [[si], count, 0]
        else:
            # update start indexes
            start_indexes = f[0]
            start_indexes.append(si)

            count = f[1]
            fst_start = start_indexes[f[2]]
            # if the word's first start index is inside my current window, increment count ((endindex is exclusive))
            if fst_start >= i - l:
                #in window
                count += 1
                frequency[word][1] = count
            else:
                #not in window, do not increment counting, update first start index (remove first)
                frequency[word][1] = count
                frequency[word][2] += 1
        if count == t:
            #if counter in this window reaches minimun, it's a winner
            clump_patterns.add(word)
        i += 1
    return clump_patterns

# print find_clump_patterns("CCGACAGGCTAGTCTATAATCCTGAGGCGTTACCCCAATACCGTTTACCGTGGGATTTGCTACTACAACTCCTGAGCGCTACATGTACGAAACCATGTTATGTAT", 4, 30, 3)
print pattern_count("TACGCATTACAAAGCACA", "AA")
# print ' '.join(str(x) for x in skew_g_c("GAGCCACCGCGATA"))
# print ' '.join(str(x) for x in minimum_skew(open("data/skew.txt").read()))
# print ' '.join(str(x) for x in minimum_skew("CATTCCAGTACTTCATGATGGCGTGAAGA"))

# print find_clump_patterns("AABABABBABA", 3, 9, 3)
# print find_clump_patterns("AGAGAGAGAGA", 3, 9, 3)
# print find_clump_patterns("AGAG", 2, 4, 2)
# print find_clump_patterns("AGAG", 1, 4, 2)
# print find_clump_patterns("AGAG", 2, 3, 2)
# print find_clump_patterns("BAG", 3, 3, 1)
# print len(find_clump_patterns(open("data/test_clump.txt").read(), 12, 595, 19))
# print ' '.join(find_clump_patterns(open("data/test_clump.txt").read(), 12, 595, 19))
# AGAGTGATTGCG GTGGATAGCCTA GTGATCCACCGA GATAGTTGGTCT ACTTCCAAACAG TACTCCTGAAGT TTGCAAACTGAC CCGCACGAAGTA ATAACGATTTCC

# print ' '.join(find_clump_patterns(open("data/test_clump2.txt").read(), 11, 566, 18))
# AAACCAGGTGG
# print ' '.join(find_clump_patterns(open("data/test_clump3.txt").read(), 9, 594, 18))

# print find_clump_patterns("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", 5, 50, 4)
# print len(find_clump_patterns(open("data/E-coli.txt").read(), 9, 500, 3))
# print len(find_clump_patterns(open("data/E-coli.txt").read(), 10, 200, 4))
# print ' '.join(find_clump_patterns(open("data/E-coli.txt").read(), 10, 200, 4))
# print len(find_clump_patterns(open("data/E-coli.txt").read(), 9, 500, 4))
# print len(find_clump_patterns(open("data/E-coli.txt").read(), 9, 500, 5))
# print len(find_clump_patterns(open("data/E-coli.txt").read(), 9, 500, 6))
# For t=4 I get 588;for t=5 I get 288;t=6 I get 34;
# print len(open("data/E-coli.txt").read()) / 500

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
