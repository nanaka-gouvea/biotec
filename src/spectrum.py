import re
from itertools import cycle

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
    # and theoretical_spectrum(p) == spec
    p_spec_size = (count_fragment_cyclic_peptide(len(p)) + 2) if cyclic else count_fragment_linear_peptide(len(p))
    if p_spec_size == len(spec) and theoretical_spectrum(p, cyclic) == spec:
        if total_mass(p) == spec[-1]:
            return [p]
        else:
            return None
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

# s_spec = re.compile("\s+").split("0	97	99	113	114	128	128	147	147	163	186	227	241	242	244 260	261	262	283	291	333	340	357	388	389	390	390	405	430	430 447	485	487	503	504	518	543	544	552	575	577	584	631	632	650 651	671	672	690	691	738	745	747	770	778	779	804	818	819	835 837	875	892	892	917	932	932	933	934	965	982	989	1031	1039	1060 1061	1062	1078	1080	1081	1095	1136	1159	1175	1175	1194	1194	1208	1209	1223 1225	1322")
# spectrum = [int(x) for x in s_spec]
# print find_cyclic_peptides_by_spectrum(spectrum)