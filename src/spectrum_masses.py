import re
from itertools import cycle
from src.tools.OrderedSet import OrderedSet

__author__ = 'natalia'

from translation import get_file
from collections import defaultdict


def a_masses_map():
    am ={line[0]: int(line[2:]) for line in get_file("/resources/integer_mass_table.txt").read().splitlines()}
    return am

a_masses = a_masses_map()

masses_set = set(a_masses.values())


def reverse_a_map():
    reverse = defaultdict(lambda: [])
    for k, v in a_masses.iteritems():
        reverse[v].append(k)
    return reverse

r_a_map = reverse_a_map()


def masses_to_pep(mpeps):
    p = ""
    for m in mpeps:
        p += r_a_map[m][0]
    return p


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
    if peptide == '':
        return 0
    return sum([int(p) for p in peptide])


def theoretical_spectrum(pep, cyclic):
    masses = [0]
    subs = fragment_cyclic_peptide(pep) if cyclic else fragment_linear_peptide(pep)
    for sub in subs:
        masses.append(total_mass(sub))
    return sorted(masses)


def peptide_score(peptide, spec, cyclic):
    aux_spec = list(spec)
    score = 0
    pts = theoretical_spectrum(peptide, cyclic)
    lp = len(pts)
    i = 0
    while i < lp:
        sp = pts[i]
        if sp in aux_spec:
            score += 1
            aux_spec.remove(sp)
            pts.remove(sp)
            lp -= 1
        else:
            i += 1
    return score


def expand_peptides_board(peps, spec, aminoacids):
    expanded = {}
    for p in peps:
        for k in aminoacids:
            new_p = '-'.join([p,str(k)])
            if p == "":
                new_p = str(k)
            #always keep the linear score, for trim
            new_p_s = new_p.split("-")
            expanded[new_p] = peptide_score(new_p_s, spec, False)
    return expanded


def old_trim(leader_board, t, ties=True):
    llen = len(leader_board)
    if llen > t:
        sorted_scores = sorted(leader_board.values())
        sorted_scores.reverse()
        # limit = sorted_scores[t]
        # if ties:
        #     for i in range(len(sorted_scores))[t:]:
        #         tie_cut = sorted_scores[i]
        #         if tie_cut < limit:
        #             sorted_scores = sorted_scores[:i]
        #             break
        # else:
        sorted_scores = sorted_scores[:t]
        aux_leader_board = leader_board.copy()
        for k,v in aux_leader_board.iteritems():
            if v in sorted_scores:
                continue
            else:
                del leader_board[k]
    return leader_board


def trim(leader_board, t, ties=True):
    llen = len(leader_board)
    if llen > t:
        sorted_peptides = sorted(leader_board.iterkeys(), key=lambda k: leader_board[k])
        sorted_peptides.reverse()
        if ties:
            limit = leader_board[sorted_peptides[t]]
            aux_leader_board = leader_board.copy()
            for k,v in aux_leader_board.iteritems():
                if v >= limit:
                    continue
                else:
                    del leader_board[k]
        else:
            leader_board = {k:leader_board[k] for k in sorted_peptides[:t]}
    return leader_board


def leader_board_experimental_spectrum(spec, t, c_masses):
    parent_mass = spec[-1]
    leader = ("",0)
    leader_board = {"":0}
    while len(leader_board) > 0:
        leader_board = expand_peptides_board(leader_board, spec, c_masses)
        for pep in leader_board.copy().keys():
            pep_s = pep.split("-")
            p_mass = total_mass(pep_s)
            if p_mass == parent_mass:
                s = peptide_score(pep_s, spec, True)
                if s >= leader[1]:
                    leader = (pep,s)
            elif p_mass > parent_mass:
                del leader_board[pep]
        trim(leader_board, t)
    return leader[0]


def leader_board_experimental_spectrum_linear_collection(spec, t, c_masses, h):
    parent_mass = spec[-1]
    leaders = {"":0}
    leader_board = {"":0}
    while len(leader_board) > 0:
        leader_board = expand_peptides_board(leader_board, spec, c_masses)
        for pep, s in leader_board.copy().iteritems():
            pep_s = pep.split("-")
            p_mass = total_mass(pep_s)
            if p_mass == parent_mass:
                leaders[pep] = s
            elif p_mass > parent_mass:
                del leader_board[pep]
        trim(leader_board, t)

    trim(leaders, h)
    return leaders


def convolution(spec):
    conv = []
    for ma in spec:
        for mb in spec:
            diff = ma - mb
            if diff > 0:
                conv.append(diff)
    return sorted(conv)


def convolution_map(spec):
    conv = {}
    for ma in spec:
        for mb in spec:
            diff = ma - mb
            if (diff >= 57) and (diff <= 200):
                try:
                    conv[diff] += 1
                except KeyError:
                    conv[diff] = 1
    return conv


def consistent(pep, spec):
    l = len(theoretical_spectrum(pep, False))
    if peptide_score(pep, spec, False) == l:
        return True
    else:
        return False


def convolution_cyclopeptide_sequencing(m, n, spectrum):
    conv = convolution_map(spectrum)
    trim(conv, m)
    return leader_board_experimental_spectrum(spectrum, n, conv.keys())


def convolution_cyclopeptide_sequencing_collection(m, n, spectrum, h):
    conv = convolution_map(spectrum)
    trim(conv, m)
    return leader_board_experimental_spectrum_linear_collection(spectrum, n, conv.keys(), h)
