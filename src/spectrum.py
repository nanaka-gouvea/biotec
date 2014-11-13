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


def expand_peptides_board(peps, spec, aminoacids):
    expanded = {}
    for p in peps:
        for k in aminoacids:
            new_p = p + k
            expanded[new_p] = peptide_score(new_p, spec, False)
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
                sorted_scores = sorted_scores[:t]
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
        leader_board = expand_peptides_board(leader_board, spec, aminoacids)
        for pep,s in leader_board.copy().iteritems():
            if total_mass(pep) == parent_mass:
                if s > leader[1]:
                    leader = (pep,s)
            elif total_mass(pep) > parent_mass:
                del leader_board[pep]
        # print len(leader_board)
        trim_peptides(leader_board, t)
        # print len(leader_board)
        # print leader_board
    return [leader[0]]


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


# spe = [int(x) for x in "0 71 71 71 87 97 97 99 101 103 113 113 114 115 128 128 129 137 147 163 163 170 184 184 186 186 190 211 215 226 226 229 231 238 241 244 246 257 257 276 277 278 299 300 312 316 317 318 318 323 328 340 343 344 347 349 356 366 370 373 374 391 401 414 414 415 419 427 427 431 437 441 446 453 462 462 462 470 472 502 503 503 511 515 529 530 533 533 540 543 547 556 559 569 574 575 584 590 600 600 604 612 616 617 630 640 640 643 646 648 660 671 683 684 687 693 703 703 719 719 719 729 730 731 737 740 741 745 747 754 774 780 784 790 797 800 806 818 826 827 832 833 838 846 846 847 850 868 869 877 884 889 893 897 903 908 913 917 930 940 947 956 960 960 961 964 965 966 983 983 985 1002 1009 1010 1011 1021 1031 1031 1036 1053 1054 1058 1059 1062 1063 1074 1076 1084 1092 1103 1113 1122 1124 1130 1133 1134 1145 1146 1146 1149 1150 1155 1156 1171 1173 1174 1187 1191 1193 1200 1212 1221 1233 1240 1242 1246 1259 1260 1262 1277 1278 1283 1284 1287 1287 1288 1299 1300 1303 1309 1311 1320 1330 1341 1349 1357 1359 1370 1371 1374 1375 1379 1380 1397 1402 1402 1412 1422 1423 1424 1431 1448 1450 1450 1467 1468 1469 1472 1473 1473 1477 1486 1493 1503 1516 1520 1525 1530 1536 1540 1544 1549 1556 1564 1565 1583 1586 1587 1587 1595 1600 1601 1606 1607 1615 1627 1633 1636 1643 1649 1653 1659 1679 1686 1688 1692 1693 1696 1702 1703 1704 1714 1714 1714 1730 1730 1740 1746 1749 1750 1762 1773 1785 1787 1790 1793 1793 1803 1816 1817 1821 1829 1833 1833 1843 1849 1858 1859 1864 1877 1886 1890 1893 1900 1900 1903 1904 1918 1922 1930 1930 1931 1961 1963 1971 1971 1971 1980 1987 1992 1996 2002 2006 2006 2014 2018 2019 2019 2032 2042 2059 2060 2063 2067 2077 2084 2086 2089 2090 2093 2105 2110 2115 2115 2116 2117 2121 2133 2134 2155 2156 2157 2176 2176 2187 2189 2192 2195 2202 2204 2207 2207 2218 2222 2243 2247 2247 2249 2249 2263 2270 2270 2286 2296 2304 2305 2305 2318 2319 2320 2320 2330 2332 2334 2336 2336 2346 2362 2362 2362 2433".split(" ")]
# print masses_leaderboard_cyclic(spe, 325)[0]

#dataset:
# print masses_to_pep([int(x) for x in "97-129-97-147-99-71-186-71-113-163-115-71-113-128-103-87-128-101-137-163-114".split("-")])
#mine:
# print masses_to_pep([int(x) for x in "71-99-147-97-129-97-114-163-137-101-128-87-103-71-57-113-71-115-163-113-71-129-57".split("-")])

# PEPFVAWAIYDAIKCSKTHYN
# AVFPEPNYHTKSCAGIADYIAEG

# print total_mass("PEPFVAWAIYDAIKCSKTHYN"), peptide_score("PEPFVAWAIYDAIKCSKTHYN", spe, False), len("PEPFVAWAIYDAIKCSKTHYN")
# print total_mass("AVFPEPNYHTKSCAGIADYIAEG"), peptide_score("AVFPEPNYHTKSCAGIADYIAEG", spe, False), len("AVFPEPNYHTKSCAGIADYIAEG")

# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 714, 2280, 2126, 1061, 439, 118, 66, 0


#this takes almost 3 minutes to run: (and gives right answer)
# spe = [int(x) for x in "0 71 87 99 99 99 99 99 103 113 114 128 128 129 129 137 147 156 163 186 186 186 198 198 215 227 228 231 234 236 241 243 246 255 257 266 269 273 276 285 300 315 318 323 326 337 340 342 345 354 356 365 368 375 384 394 397 401 410 420 422 429 429 439 455 462 465 467 473 474 479 481 496 496 504 509 519 523 528 538 540 551 552 558 561 576 586 595 595 595 602 611 618 622 637 638 641 651 653 657 660 665 665 667 675 694 694 705 715 721 724 738 738 740 742 750 752 766 774 774 781 788 793 794 804 804 804 837 841 849 851 852 869 871 875 877 880 880 887 901 902 903 903 903 924 936 938 940 941 951 970 974 979 980 983 990 1002 1005 1015 1015 1023 1031 1032 1040 1050 1057 1061 1066 1067 1069 1092 1103 1114 1118 1118 1122 1126 1127 1130 1139 1144 1146 1156 1160 1160 1169 1171 1195 1197 1205 1213 1214 1217 1217 1226 1243 1246 1255 1255 1259 1259 1267 1270 1278 1289 1298 1300 1313 1316 1325 1326 1332 1333 1342 1345 1358 1360 1369 1380 1388 1391 1399 1399 1403 1403 1412 1415 1432 1441 1441 1444 1445 1453 1461 1463 1487 1489 1498 1498 1502 1512 1514 1519 1528 1531 1532 1536 1540 1540 1544 1555 1566 1589 1591 1592 1597 1601 1608 1618 1626 1627 1635 1643 1643 1653 1656 1668 1675 1678 1679 1684 1688 1707 1717 1718 1720 1722 1734 1755 1755 1755 1756 1757 1771 1778 1778 1781 1783 1787 1789 1806 1807 1809 1817 1821 1854 1854 1854 1864 1865 1870 1877 1884 1884 1892 1906 1908 1916 1918 1920 1920 1934 1937 1943 1953 1964 1964 1983 1991 1993 1993 1998 2001 2005 2007 2017 2020 2021 2036 2040 2047 2056 2063 2063 2063 2072 2082 2097 2100 2106 2107 2118 2120 2130 2135 2139 2149 2154 2162 2162 2177 2179 2184 2185 2191 2193 2203 2219 2229 2229 2236 2238 2248 2257 2261 2264 2274 2283 2290 2293 2302 2304 2313 2316 2318 2321 2332 2335 2340 2343 2358 2373 2382 2385 2389 2392 2401 2403 2412 2415 2417 2422 2424 2427 2430 2431 2443 2460 2460 2472 2472 2472 2495 2502 2511 2521 2529 2529 2530 2530 2544 2545 2555 2559 2559 2559 2559 2559 2571 2587 2658".split(" ")]
# print masses_leaderboard_cyclic(spe, 365)[0]
# 99-99-156-113-128-99-99-147-129-186-57-57-129-99-137-186-87-128-103-163-71-99-87