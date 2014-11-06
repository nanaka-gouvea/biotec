import os
from collections import defaultdict
from frequent_words import *

__author__ = 'natalia'


def get_file(filepath):
    return open(os.path.dirname(__file__) + filepath)


def t_map():
    return {line[0:3]:line[4:] if len(line) > 3 else " " for line in get_file("/resources/genetic_code.txt").read().splitlines()}


def amino_map():
    return {line[0]:line[2:] if len(line) > 2 else "" for line in get_file("/resources/aminoacids.txt").read().splitlines()}

t_map = t_map()


def reverse_t_map():
    reverse = defaultdict(lambda: [])
    for k, v in t_map.iteritems():
        reverse[v].append(k)
    return reverse

reverse_t_map = reverse_t_map()


def translate_peptide(seq):
    aminoacids = [t_map[k] for k in [seq[i:i + 3] for i in range(0, len(seq), 3)]]
    return ''.join(aminoacids)


def count_linear_peptyde_codons(aminoacids):
    count = 1
    for a in aminoacids:
        count *= len(reverse_t_map[a])
    return count


def transcribe(genome):
    return genome.replace("T", "U")


def unscribe(genome):
    return genome.replace("U", "T")


def verify_encoding(genome, i, k, peptide, seqs):
    while i <= len(genome) - k:
        word = genome[i:i + k]
        transc1 = transcribe(word)
        transc2 = transcribe(reverse_complement(word))
        if translate_peptide(transc1) == peptide or translate_peptide(transc2) == peptide:
            seqs.append(word)
        i += 1


def find_peptide_encoding(genome, peptide):
    k = len(peptide) * 3
    seqs = []
    for i in [0,1,2]:
        verify_encoding(genome, i, k, peptide, seqs)
    return seqs


def better_verify_encoding(genome, k, peptide):
    subs = []
    for i in [0, 1, 2]:
        e = -((len(genome) - i) % 3)
        seq = genome[i:e if e != 0 else None]
        # todo better performance, add matches while translating, stop translating (and jump to next k)at each aminoacid mismatch
        translation = translate_peptide(seq)
        for j in range(len(translation) - k + 1):
            word = translation[j:j + k]
            if word == peptide:
                subs.append(unscribe(seq[j*3:(j*3)+k*3]))
    return subs


def better_find_peptide_encoding(genome, peptide):
    k = len(peptide)
    subs = []
    subs.extend(better_verify_encoding(transcribe(genome), k, peptide))
    subs.extend([reverse_complement(s) for s in better_verify_encoding(transcribe(reverse_complement(genome)), k, peptide)])
    return subs

# gen = ''.join(get_file("/data/B_brevis.txt").read().splitlines())
gen = get_file("/data/B_brevis.txt").read().replace("\n", "")
tyrocidine = get_file("/resources/Tyrocidine_B1.txt").readline().strip()
# print len(better_find_peptide_encoding(gen, tyrocidine))