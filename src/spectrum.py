import re
from itertools import cycle
from src.tools.OrderedSet import OrderedSet

__author__ = 'natalia'

from translation import get_file
from collections import defaultdict


def a_masses_map():
    am ={line[0]: int(line[2:]) for line in get_file("/resources/integer_mass_table.txt").read().splitlines()}
    # am['']=0
    return am

a_masses = a_masses_map()
# a_masses_0 = a_masses_map().copy()
# a_masses_0['']=0


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


def masses_find_cyclic_peptides_by_spectrum(spec):
    u_m_a = list(a_masses.keys())
    #removes aminoacids which have homonimal masses
    u_m_a.remove("L")
    u_m_a.remove("Q")
    possible = possible_peptides_from_spectrum(spec, "", True, u_m_a)
    result = set()
    for p in possible:
        result.add("-".join([str(s) for s in [a_masses[a] for a in p]]))
    return result

#to linear cyclic should be false
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


def expand_peptides_board(peps, spec, aminoacids, p_mass):
    expanded = {}
    for p in peps:
        for k in aminoacids:
            new_p = p + k
            n_mass = total_mass(new_p)
            if n_mass == p_mass:
                expanded[new_p] = peptide_score(new_p, spec, True)
            elif n_mass < p_mass:
                expanded[new_p] = peptide_score(new_p, spec, False)
                #does not expand to peptides with mass greater than the parent mass
    return expanded


def trim_peptides(leader_board, t):
    llen = len(leader_board)
    if llen > t:
        sorted_scores = sorted(leader_board.values())
        sorted_scores.reverse()
        limit = sorted_scores[t]
        for i in range(len(sorted_scores))[t:]:
            low_s = sorted_scores[i]
            if low_s < limit:
                sorted_scores = sorted_scores[:t+1]
                break
        aux_leader_board = leader_board.copy()
        for k,v in aux_leader_board.iteritems():
            if v in sorted_scores:
                continue
            else:
                del leader_board[k]
    return leader_board


def leader_board_experimental_spectrum(spec, aminoacids, t):
    "t is tolerance for the leaderboard"
    parent_mass = spec[-1]
    leader = ("",0)
    leader_board = {"":0}
    while len(leader_board) > 0:
        leader_board = expand_peptides_board(leader_board, spec, aminoacids, parent_mass)
        for pep,s in leader_board.copy().iteritems():
            if total_mass(pep) == parent_mass:
                if s > leader[1]:
                    leader = (pep,s)
                    # elif total_mass(pep) > parent_mass:
                    #     del leader_board[pep]
        # print len(leader_board)
        trim_peptides(leader_board, t)
        # print len(leader_board)
        # print leader_board
    return [leader[0]]

def leader_board_experimental_spectrum_many(spec, aminoacids, t):
    "t is tolerance for the leaderboard"
    parent_mass = spec[-1]
    leader = [("",0)]
    leader_board = {"":0}
    while len(leader_board) > 0:
        leader_board = expand_peptides_board(leader_board, spec, aminoacids)
        for pep,s in leader_board.copy().iteritems():
            if total_mass(pep) == parent_mass:
                if s > leader[0][1]:
                    leader = [(pep,s)]
                elif s == leader[0][1]:
                    leader.append((pep,s))
            elif total_mass(pep) > parent_mass:
                del leader_board[pep]
        trim_peptides(leader_board, t)
    return leader


def masses_leaderboard_cyclic(spec, n):
    u_m_a = list(a_masses.keys())
    #removes aminoacids which have homonimal masses
    u_m_a.remove("L")
    u_m_a.remove("Q")
    possible = leader_board_experimental_spectrum(spec, u_m_a, n)
    result = []
    for p in possible:
        result.append("-".join([str(s) for s in [a_masses[a] for a in p]]))
    return result


def convolution(spec):
    conv = []
    for ma in spec:
        for mb in spec:
            diff = ma - mb
            #todo for better performance, filter 57 <x>200 here
            if diff > 0:
                conv.append(diff)
    return conv

    # teste = [str(y) for y in convolution([int(x) for x in "465 473 998 257 0 385 664 707 147 929 87 450 748 938 998 768 234 722 851 113 700 957 265 284 250 137 317 801 128 820 321 612 956 434 534 621 651 129 421 337 216 699 347 101 464 601 87 563 738 635 386 972 620 851 948 200 156 571 551 522 828 984 514 378 363 484 855 869 835 234 1085 764 230 885".split(" ")])]


spe = [int(m) for m in "0 71 71 71 87 97 97 99 101 103 113 113 114 115 128 128 129 137 147 163 163 170 184 184 186 186 190 211 215 226 226 229 231 238 241 244 246 257 257 276 277 278 299 300 312 316 317 318 318 323 328 340 343 344 347 349 356 366 370 373 374 391 401 414 414 415 419 427 427 431 437 441 446 453 462 462 462 470 472 502 503 503 511 515 529 530 533 533 540 543 547 556 559 569 574 575 584 590 600 600 604 612 616 617 630 640 640 643 646 648 660 671 683 684 687 693 703 703 719 719 719 729 730 731 737 740 741 745 747 754 774 780 784 790 797 800 806 818 826 827 832 833 838 846 846 847 850 868 869 877 884 889 893 897 903 908 913 917 930 940 947 956 960 960 961 964 965 966 983 983 985 1002 1009 1010 1011 1021 1031 1031 1036 1053 1054 1058 1059 1062 1063 1074 1076 1084 1092 1103 1113 1122 1124 1130 1133 1134 1145 1146 1146 1149 1150 1155 1156 1171 1173 1174 1187 1191 1193 1200 1212 1221 1233 1240 1242 1246 1259 1260 1262 1277 1278 1283 1284 1287 1287 1288 1299 1300 1303 1309 1311 1320 1330 1341 1349 1357 1359 1370 1371 1374 1375 1379 1380 1397 1402 1402 1412 1422 1423 1424 1431 1448 1450 1450 1467 1468 1469 1472 1473 1473 1477 1486 1493 1503 1516 1520 1525 1530 1536 1540 1544 1549 1556 1564 1565 1583 1586 1587 1587 1595 1600 1601 1606 1607 1615 1627 1633 1636 1643 1649 1653 1659 1679 1686 1688 1692 1693 1696 1702 1703 1704 1714 1714 1714 1730 1730 1740 1746 1749 1750 1762 1773 1785 1787 1790 1793 1793 1803 1816 1817 1821 1829 1833 1833 1843 1849 1858 1859 1864 1877 1886 1890 1893 1900 1900 1903 1904 1918 1922 1930 1930 1931 1961 1963 1971 1971 1971 1980 1987 1992 1996 2002 2006 2006 2014 2018 2019 2019 2032 2042 2059 2060 2063 2067 2077 2084 2086 2089 2090 2093 2105 2110 2115 2115 2116 2117 2121 2133 2134 2155 2156 2157 2176 2176 2187 2189 2192 2195 2202 2204 2207 2207 2218 2222 2243 2247 2247 2249 2249 2263 2270 2270 2286 2296 2304 2305 2305 2318 2319 2320 2320 2330 2332 2334 2336 2336 2346 2362 2362 2362 2433".split(" ")]
teste = masses_leaderboard_cyclic(spe, 325)
# 97-129-97-147-99-71-186-71-113-163-115-71-113-128-103-87-128-101-137-163-114
# 97-114-163-137-101-128-87-103-128-113-71-115-163-113
# 97-129-97-114-163-137-101-128-87-103-128-113-71-115-163-113
# 71-186-71-99-147-97-129-97-114-163-137-101-128-87-103-128-113-71-115-163-113
# 71-186-71-99-147-97-129-97-114-163-137-101-128-87-103-128-113-71-115-163-113