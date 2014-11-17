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


def trim(leader_board, t):
    llen = len(leader_board)
    if llen > t:
        sorted_scores = sorted(leader_board.values())
        sorted_scores.reverse()
        limit = sorted_scores[t]
        # sorted_scores[96]
        for i in range(len(sorted_scores))[t:]:
            tie_cut = sorted_scores[i]
            if tie_cut < limit:
                sorted_scores = sorted_scores[:i]
                break
        aux_leader_board = leader_board.copy()
        for k,v in aux_leader_board.iteritems():
            if v in sorted_scores:
                continue
            else:
                del leader_board[k]
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


# spe = [int(x) for x in "0 57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493".split(" ")]
# print convolution_cyclopeptide_sequencing(20, 60, spe)

spe = sorted([int(x) for x in "1023 502 1285 99 1262 1285 986 495 1325 691 783 1186 542 1236 236 1422 1037 186 385 436 773 300 236 1031 1365 231 937 1259 486 792 1186 1323 936 186 1210 288 771 1236 163 57 979 824 283 1024 250 451 1222 1123 1073 485 99 1073 735 391 938 356 887 394 1049 536 834 200 649 651 880 1036 1139 1156 1191 349 97 674 571 851 1165 1019 634 299 835 767 196 266 299 160 488 103 398 672 1134 103 790 349 1028 598 403 386 128 655 974 920 1040 1137 1237 0 1294 535 837 585 632 788 868 448 630 1309 137 887 750 670 1123 886 971 288 1134 382 535 1066 687 588 113 484 185 587 443 1172 99 934 1335 1122 1323 554 137 752 285 87 1319 748 927 1319 731 1323 373 891 1226 639 531 257 212 399".split(" ")])
print convolution_cyclopeptide_sequencing(16, 326, spe)


