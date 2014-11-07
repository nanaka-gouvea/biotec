__author__ = 'natalia'

from translation import *


def a_masses_map():
    return {line[0]:int(line[2:]) for line in get_file("/resources/integer_mass_table.txt").read().splitlines()}

a_masses = a_masses_map()


def reverse_a_map():
    reverse = defaultdict(lambda: [])
    for k, v in a_masses.iteritems():
        reverse[v].append(k)
    return reverse

r_a_map = reverse_a_map()


def count_fragment_peptide(l):
    return l * (l - 1)


def fragment_peptide(pep):
    subs = [pep]
    c_pep = pep + pep
    for i in range(len(pep)):
        for l in range(len(pep))[1:]:
            subs.append(c_pep[i:i + l])
    return subs


def theoretical_spectrum(pep):
    masses = [0]
    for sub in fragment_peptide(pep):
        masses.append(sum([a_masses[a] for a in sub]))
    return sorted(masses)

