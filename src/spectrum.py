import re
from itertools import cycle
from src.tools.OrderedSet import OrderedSet

__author__ = 'natalia'

from translation import get_file
from collections import defaultdict


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
    # adds one for the 'empty' subpeptyde
    return sum(range(l + 1)[1:]) + 1


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


def theoretical_spectrum(pep, cyclic):
    masses = [0]
    subs = fragment_cyclic_peptide(pep) if cyclic else fragment_linear_peptide(pep)
    for sub in subs:
        masses.append(total_mass(sub))
    return sorted(masses)


def possible_peptides_from_spectrum(spec, p, cyclic, aminoacids=None):
    p_spec_size = (count_fragment_cyclic_peptide(len(p)) + 2) if cyclic else count_fragment_linear_peptide(len(p))
    if p_spec_size == len(spec):
        if theoretical_spectrum(p, cyclic) == spec:
            if total_mass(p) == spec[-1]:
                return [p]
            else:
                return []
        else:
            return []
    branch = []
    if aminoacids is None:
        aminoacids = a_masses.keys()
    b_aminoacids = list(aminoacids)
    for k in aminoacids:
        new = p + k
        if total_mass(new) in spec:
            branch.extend(possible_peptides_from_spectrum(spec, new, cyclic, b_aminoacids))
        else:
            if len(new) == 1:
                b_aminoacids.remove(k)
    return branch


def find_cyclic_peptides_by_spectrum(spec):
    return possible_peptides_from_spectrum(spec, "", True)


def find_linear_peptides_by_spectrum(spec):
    return possible_peptides_from_spectrum(spec, "", False)


def masses_possible_peptides_from_spectrum(spec, p, cyclic, amino_masses=None):
    p_spec_size = (count_fragment_cyclic_peptide(len(p)) + 2) if cyclic else count_fragment_linear_peptide(len(p))
    if p_spec_size == len(spec):
        if theoretical_spectrum(p, cyclic) == spec:
            if sum(p) == spec[-1]:
                return p
            else:
                return []
        else:
            return []
    branch = []
    if amino_masses is None:
        amino_masses = [str(x) for x in set(a_masses.values())]
    b_aminoacids = list(amino_masses)
    for k in amino_masses:
        new = p.append(k)
        if sum(new) in spec:
            branch.extend(possible_peptides_from_spectrum(spec, new, cyclic, b_aminoacids))
        else:
            if len(new) == 1:
                b_aminoacids.remove(k)
    return branch


def masses_find_cyclic_peptides_by_spectrum(spec):
    possible = possible_peptides_from_spectrum(spec, "", True)
    result = set()
    for p in possible:
        result.add("-".join([str(s) for s in [a_masses[a] for a in p]]))
    return result


def peptide_score(peptide, spec, cyclic):
    score = 0
    pts = theoretical_spectrum(peptide, cyclic)
    lp = len(pts)
    i = 0
    while i < lp:
        sp = pts[i]
        if sp in spec:
            score += 1
            #TODO melhorar performance
            spec.remove(sp)
            pts.remove(sp)
            lp -= 1
        else:
            i += 1
    return score


def masses_possible_peptides_from_experimental_spectrum(spec, p, cyclic, amino_masses=None):
    p_spec_size = (count_fragment_cyclic_peptide(len(p)) + 2) if cyclic else count_fragment_linear_peptide(len(p))
    if p_spec_size == len(spec):
        if peptide_score(p, spec, cyclic) >= 1:
            if sum(p) == spec[-1]:
                return p
            else:
                return []
        else:
            return []
    branch = []
    if amino_masses is None:
        amino_masses = [str(x) for x in set(a_masses.values())]
    b_aminoacids = list(amino_masses)
    for k in amino_masses:
        new = p.append(k)
        if sum(new) in spec:
            branch.extend(possible_peptides_from_spectrum(spec, new, cyclic, b_aminoacids))
        else:
            if len(new) == 1:
                b_aminoacids.remove(k)
    return branch

