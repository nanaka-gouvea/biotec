import os
from collections import defaultdict
from frequent_words import *

__author__ = 'natalia'

def get_file(filepath):
    return open(os.path.dirname(__file__) + filepath)


def t_map():
    return {line[0:3]:line[4:-1] if len(line) > 3 else "" for line in get_file("/resources/genetic_code.txt").readlines()}


def amino_map():
    return {line[0]:line[2:-1] if len(line) > 2 else "" for line in get_file("/resources/aminoacids.txt").readlines()}

amino_map = amino_map()
t_map = t_map()


def reverse_t_map():
    reverse = defaultdict(lambda: [])
    for k, v in t_map.iteritems():
        reverse[v].append(k)
    return reverse

reverse_t_map = reverse_t_map()


def translate_peptide(seq):
    peptides = [t_map[k] for k in [seq[i:i + 3] for i in range(0, len(seq), 3)]]
    return ''.join(peptides)


def count_linear_peptyde_codons(aminoacids):
    count = 1
    for a in aminoacids:
        count *= len(reverse_t_map[a])
    return count

print count_linear_peptyde_codons("KVLFPWFNQY")

print reverse_complement("ATGC")
