# coding=utf-8
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
    k = len(pattern)
    # end index is exclusive
    for i in range(len(text) - k + 1):
        if text[i:i + k] == pattern:
            starts.append(i)
    return starts

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


def hamming_d(a, b):
    i = 0
    count = 0
    while i < len(a):
        if a[i] != b[i]:
            count += 1
        i += 1
    return count


def pattern_starts_aprox(pattern, text, d):
    starts = []
    for i in range(len(text) - (len(pattern) - 1)):
        word = text[i:i + len(pattern)]
        if hamming_d(word, pattern) <= d:
            starts.append(i)
    return ' '.join(str(x) for x in starts)


def pattern_count_aprox(text, pattern, d):
    count = 0
    for i in range(len(text) - (len(pattern) - 1)):
        word = text[i:i + len(pattern)]
        if hamming_d(word, pattern) <= d:
            count += 1
    return count


def old_most_frequent_words_aprox(text, size, d):
    frequency = {}
    max_count = 0
    frequents = []
    for i in range(len(text) - (size - 1)):
        word = text[i:i + size]
        if word not in frequency:
            frequency[word] = 0
        matches = False
        for w in frequency.keys():
            if hamming_d(w, word) <= d:
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
    return result


def arrange_repeated_dna(size):
    return arrange_repeated(pool, size)


def generate_neighbourhood(letters, size, pattern, d):
    if size == 0:
        return []
    result = []
    if size == 1:
        for l in letters:
            result.append(l)
        return result
    for a in generate_neighbourhood(letters, size - 1, pattern, d):
        for i in range(len(letters)):
            nei = letters[i] + a[:i] + a[i:]
            if len(pattern) != len(nei) or hamming_d(pattern, nei) <= d:
                result.append(nei)
    return result


def better_neighbors_reverse(letters, size, pattern, d):
    if size == 0:
        return []
    result = []
    if size == 1:
        for l in letters:
            result.append(l)
        return result
    for a in generate_neighbourhood(letters, size - 1, pattern, d):
        for i in range(len(letters)):
            nei = letters[i] + a[:i] + a[i:]
            if len(pattern) != len(nei) or hamming_d(pattern, nei) <= d:
                result.append(nei)
                result.append(reverse_complement(nei))
    return result

pool = "ACGT"
def pattern_to_number(pattern):
    possible = arrange_repeated(pool, len(pattern))
    for i in range(len(possible)):
        if possible[i] == pattern:
            return i


def number_to_pattern(i, size):
    possible = arrange_repeated(pool, size)
    return possible[i]


def frequency_array(seq, size):
    possible = arrange_repeated(pool, size)
    fa = collections.OrderedDict()
    for p in possible:
        fa[p] = 0
    for i in range(len(seq) - (size - 1)):
        fa[seq[i:i + size]] += 1
    return fa.values()


def clump_finding(genome, k, l, t):
    start_indexes = {}
    clump_patterns = set()
    i = k
    while i <= len(genome):
        current_start_ix = i - k
        word = genome[current_start_ix:i]
        try:
            start_indexes[word].append(current_start_ix)
        except KeyError:
            start_indexes[word] = [current_start_ix]
        word_starting_indexes = start_indexes[word]
        len_ixs = len(word_starting_indexes)
        if len_ixs >= t:
            #start index of the first match in the last t matches
            last_tth_start_ix = word_starting_indexes[len_ixs - t]
            if i - last_tth_start_ix <= l:
                clump_patterns.add(word)
        i += 1
    return clump_patterns


def smart_neighbourhood_count(seq, max_hamming):
    size = len(seq)
    return pow(len(pool),size) - pow(len(pool), size - max_hamming)

def neighbourhood_count(seq, max_hamming):
    count = 0
    for i in arrange_repeated(pool, len(seq)):
        if hamming_d(i, seq) <= max_hamming:
            count += 1
    return count

def neighbors(pattern, max_hamming, possible=None):
    neighbours = []
    if possible is None:
        possible = arrange_repeated(pool, len(pattern))
    for i in possible:
        if hamming_d(i, pattern) <= max_hamming:
            neighbours.append(i)
    return neighbours

# print smart_neighbourhood_count("ACG", 1)
# print neighbourhood_count("ACG", 1)

def find_specific_clump():
    sixs = pattern_starts("CGCATCCGG", open("data/E-coli.txt").read())
    print sixs
    i = 4
    while i < (len(sixs) - 4):
        if (sixs[i + 4] - sixs[i]) <= 601:
            print [sixs[i], sixs[i + 1], sixs[i + 2], sixs[i + 3], sixs[i + 4]]
        i += 1


def most_frequent_words_aprox(seq, k, d):
    most = ()
    n_hood = []

    for i in range(len(seq) - (k - 1)):
        n_hood.extend(generate_neighbourhood(pool, k, seq[i:i + k], d))

    n_hood = sorted(n_hood)
    max_count = 0
    count = 1
    for i in range(len(n_hood) - 1):
        word = n_hood[i]
        if word == n_hood[i + 1]:
            count += 1
            if count > max_count:
                max_count = count
                most = [word]
            elif count == max_count:
                most.append(word)
        else:
            count = 1
    return most


def most_frequent_words_aprox_reverse(seq, k, d):
    most = ()
    n_hood = []

    for i in range(len(seq) - (k - 1)):
        pattern = seq[i:i + k]
        n_hood.extend(better_neighbors_reverse(pool, k, pattern, d))

    n_hood = sorted(n_hood)
    max_count = 0
    count = 1
    for i in range(len(n_hood) - 1):
        word = n_hood[i]
        if word == n_hood[i + 1]:
            count += 1
            if count > max_count:
                max_count = count
                most = [word]
            elif count == max_count:
                most.append(word)
        else:
            count = 1
    return most