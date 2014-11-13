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


def expand_peptides_board(peps, spec):
    expanded = {}
    for p in peps:
        for k in a_masses.keys():
            expanded[p+k] = peptide_score(p, spec, False)
    return expanded


# Trim(Leaderboard, Spectrum, N, AminoAcid, AminoAcidMass)
# for j = 1 to |Leaderboard
#         Peptide = j-th peptide in Leaderboard
# LinearScores(j) = LinearScore(Peptide, Spectrum)
# sort Leaderboard according to the decreasing order of scores in LinearScores
# sort LinearScores in decreasing order
# for j = N + 1 to |Leaderboard|
#   if LinearScores(j) < LinearScores(N)
#     remove all peptides starting from the j-th peptide from Leaderboard
# return Leaderboard


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

#     LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N)
#     Leaderboard = {empty peptide}
#     LeaderPeptide = empty peptide
#     while Leaderboard is non-empty
#         Leaderboard = Expand(Leaderboard)
#         for each Peptide in Leaderboard
#             if Mass(Peptide) = ParentMass(Spectrum)
#               if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
#                 LeaderPeptide = Peptide
#             else if Mass(Peptide) > ParentMass(Spectrum)
#               remove Peptide from Leaderboard
#         Leaderboard = Trim(Leaderboard, Spectrum, N)
#     output LeaderPeptide

def leader_board_experimental_spectrum(spec, cyclic, aminoacids, t):
    "t is tolerance for the leaderboard"
    parent_mass = spec[-1]
    leader = ""
    leader_board = {leader:0}
    while len(leader_board) > 0:
        # leader_board = expand_peptides(leader_board)
        for pep in leader_board:
            if total_mass(pep) == parent_mass:
                #if
                return leader
    # if total_mass(p) == spec[-1]:
    #     return [p]
    # branch = [""]
    # for b in branch:
    #     n_branches = {}
    #     for k in aminoacids:
    #         new = b + k
    #         n_branches[new] = linear_peptide_score(new, spec, cyclic)
    #     branch.extend(trim_peptides(n_branches, t))
    return leader


def masses_find_cyclic_peptides_by_experimental_spectrum(spec):
    u_m_a = list(a_masses.keys())
    #removes aminoacids which have homonimal masses
    u_m_a.remove("L")
    u_m_a.remove("Q")
    possible = leader_board_experimental_spectrum(spec, True, u_m_a, 7)
    result = set()
    for p in possible:
        result.add("-".join([str(s) for s in [a_masses[a] for a in p]]))
    return result


