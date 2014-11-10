__author__ = 'natalia'

from translation import *


def a_masses_map():
    return {line[0]: int(line[2:]) for line in get_file("/resources/integer_mass_table.txt").read().splitlines()}

a_masses = a_masses_map()


def reverse_a_map():
    reverse = defaultdict(lambda: [])
    for k, v in a_masses.iteritems():
        reverse[v].append(k)
    return reverse

r_a_map = reverse_a_map()


def count_fragment_linear_peptide(l):
    return sum(range(l + 1)[1:])


def fragment_linear_peptide(pep):
    subs = []
    for i in range(len(pep)):
        for l in range(len(pep) - i + 1)[1:]:
            subs.append(pep[i:i + l])
    return subs


def count_fragment_cyclic_peptide(l):
    return l * (l - 1)


def fragment_cyclic_peptide(pep):
    subs = [pep]
    c_pep = pep + pep
    for i in range(len(pep)):
        for l in range(len(pep))[1:]:
            subs.append(c_pep[i:i + l])
    return subs


def total_mass(peptide):
    return sum([a_masses[a] for a in peptide])


def theoretical_spectrum(pep):
    masses = [0]
    for sub in fragment_cyclic_peptide(pep):
        masses.append(total_mass(sub))
    return sorted(masses)


def possible_peptides_from_spectrum(spec, p):
    # and theoretical_spectrum(p) == spec
    p_spec_size = count_fragment_cyclic_peptide(len(p)) + 2
    if p_spec_size == len(spec) and theoretical_spectrum(p) == spec:
        if total_mass(p) == spec[-1]:
            return [p]
        else:
            return None
    branch = []
    for k in a_masses.keys():
        new = p + k
        if total_mass(new) in spec:
            branch.extend(possible_peptides_from_spectrum(spec, new))
    return branch


def find_peptides_by_spectrum(spec):
    return possible_peptides_from_spectrum(spec, "")

# print find_peptides_by_spectrum([0, 57, 71, 128])

print fragment_linear_peptide("NQEL")