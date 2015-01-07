from sys import maxint
from translation import get_file
from translation import get_file_w
from genome import non_branching_paths
import re
from random import choice
from copy import deepcopy
from origin import hamming_d

__author__ = 'natalia'

#TODO go to pyhton 3.4 to use enums!! and then put all into classes =(
#TODO GET RID OF THOSE REPEATED OPTIMAL PATHS!!!
#todo the path choice is all f* up. there are many more choices, you cant delete the extras before all uses if it

def construct_alignment_all(l, c, trace, s1, s2, op, choice_track, retrack):
    should_retrack = False
    if len(choice_track) > 0:
        path = choice_track[-1][0]
        current = choice_track[-1][1]
        retrack.append((current, trace[(current, 0)][0]))
        del trace[(current, 0)][0]
        if len(trace[current]) == 1:
            del choice_track[-1]
            should_retrack = True
    else:
        current = ((l, c), 0)
        path = ["", ""]
    while current[0] != (0, 0):
        all_previous = trace[current]
        previous = all_previous[0]
        direction = all_previous[0][1]
        if len(all_previous) > 1:
            op += 1
            choice_track.append([path[:], current])
        #middle matrix is diagonal
        if direction == 0:
            path[0] += s1[current[0][1] - 1]
            path[1] += s2[current[0][0] - 1]
        #lower matrix is down
        elif direction == -1:
            path[0] += "-"
            path[1] += s2[current[0][0] - 1]
        #upper matrix is right
        elif direction == 1:
            path[0] += s1[current[0][1] - 1]
            path[1] += "-"
        else:
            #free taxi ride
            pass
        current = previous
    if should_retrack:
        for r in retrack[:]:
            trace.setdefault(r[0], []).append(r[1])
            retrack.remove(r)
    distance = hamming_d(path[0], path[1])
    return [path[0][::-1], path[1][::-1], distance], op


def find_all_alignments(column, line, trace, s1, s2):
    paths = []
    options = 1
    choice_track = []
    retrack = []
    # result = construct_alignment_all(line, column, trace, s1, s2, options, choice_track, retrack)
    # paths.append(result[0])
    while options > 0:
        result = construct_alignment_all(line, column, trace, s1, s2, options, choice_track, retrack)
        paths.append(result[0])
        options = result[1]
        options -= 1
    return paths


def sequence_alignment_scored(s1, s2, p, ep, subm, atype="global"):
    """
    :param s1: sequence 1
    :param s2: sequence 2
    :param p: gap opening penalty
    :param ep: gap extension penalty
    :param subm substution matrix for scoring
    :return: possible alignments with greatest score
    """
    local = atype == "local"
    globalt = atype == "global"
    fitting = atype == "fitting"
    overlap = atype == "overlap"

    line = len(s2)
    column = len(s1)
    range_c = range(column + 1)
    range_l = range(line + 1)
    #trace example: {(1,2):[((0,2),"down"),((1,1),"right")]}
    #trace example: {[1,2,-1]:[((0,2),-1),((1,1),1)]}
    #this means that the max value in node 1,2, on the lower matrix, can come from 0,2 if going down, n' from 1,1 going right
    trace = {}
    middle = []
    upper = []
    lower = []

    if globalt:
        lower = [[0]]
        lower += [[-p + (l * (-ep))] for l in range_l[:-1]]
        lower[0] += ([-p + (c * (-ep)) for c in range_c[:-1]])

        middle = [[0]]
        middle += [[-p + (l * (-ep))] for l in range_l[:-1]]
        middle[0] += ([-p + (c * (-ep)) for c in range_c[:-1]])

        upper = [[0]]
        upper += [[-p + (l * (-ep))] for l in range_l[:-1]]
        upper[0] += ([-p + (c * (-ep)) for c in range_c[:-1]])

    if globalt:
        #column 0
        for i in range_l[1:]:
            trace.setdefault(((i,0), -1), []).append(((i - 1, 0), -1))
        #line 0
        for i in range_c[1:]:
            trace.setdefault(((0,i), 1), []).append(((0, i - 1), 1))

    print "start matrix: "
    for m in middle:
        print m

    #for keeping the greatest value (sink) to local alignment
    end_line = 0
    end_column = 0

    for i in range(len(lower))[1:]:
        for j in range(len(upper[0]))[1:]:
            le = lower[i - 1][j] - ep
            lo = middle[i - 1][j] - p
            if i == 1:
                le = lo
            lowermax = max(le, lo)
            lower[i].append(lowermax)
            ue = upper[i][j - 1] - ep
            uo = middle[i][j - 1] - p
            if j == 1:
                ue = uo
            uppermax = max(ue, uo)
            upper[i].append(uppermax)

            sub_score = subm[s1[j - 1]][s2[i - 1]]
            diagonal = middle[i - 1][j - 1] + sub_score
            middlemax = max(lowermax, uppermax, diagonal)
            middle[i].append(middlemax)

            if middlemax == diagonal:
                try:
                    _ = trace[((i - 1, j - 1), 0)]
                    trace.setdefault(((i,j), 0), []).append(((i - 1, j - 1), 0))
                except KeyError:
                    try:
                        _ = trace[((i - 1, j - 1), -1)]
                        trace.setdefault(((i,j), 0), []).append(((i - 1, j - 1), -1))
                    except KeyError:
                        trace.setdefault(((i,j), 0), []).append(((i - 1, j - 1), 1))
            if middlemax == lowermax:
                if lowermax == lo:
                    trace.setdefault(((i,j), -1), []).append(((i - 1,j), 0))
                else:
                    trace.setdefault(((i,j), -1), []).append(((i - 1,j), -1))
            if middlemax == uppermax:
                if uppermax == uo:
                    trace.setdefault(((i,j), 1), []).append(((i,j - 1), 0))
                else:
                    trace.setdefault(((i,j), 1), []).append(((i,j - 1), 1))

    print "final matrix: "
    for m in middle:
        print m
    if globalt:
        end_line = line
        end_column = column
    print "FINAL SCORE: " + str(middle[end_line][end_column]) + " (" + str(end_line) + "," + str(end_column) + ")"
    return find_all_alignments(end_column, end_line, trace, s1, s2)


def read_subm(subm_name=None):
    if subm_name is None:
        mfile = get_file("/resources/pam250.txt")
    else:
        mfile = get_file("/resources/" + subm_name + ".txt")
    lines = [l.split() for l in mfile.read().splitlines()]
    xline = lines[0]
    subm = {x:{} for x in xline}
    for l in lines[1:]:
        for i in range(len(l))[1:]:
            subm[xline[i - 1]][l[0]] = int(l[i])
    return subm


def output_alignment_scored(s1, s2, p, ep, atype="global", subm_name=None):
    alignments = sequence_alignment_scored(s1, s2, p, ep, read_subm(subm_name), atype)
    print "Possible:"
    count = 0
    for a in alignments:
        count += 1
        print count, ": "
        print a[0]
        print a[1]
        print "edit distance: " + str(a[2])


# t1 = "PRTEINS"
# t2 = "PRTWPSEIN"
t1 = "YHFDVPDCWAHRYWVENPQAIAQMEQICFNWFPSMMMKQPHVFKVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE"
t2 = "YHEDVAHEDAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPIISATCARMRVRTVWE"
output_alignment_scored(t1, t2, 11, 1, "global", "blossum62")
# construct_alignment_graph(t1, t2, -5, read_subm())



#179
# YHEDV----AH----ED--AIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPI----ISATCARMRVRTVWE
# YHFDVPDCWAHRYWVENPQAIAQM-------EQICFNWFPSMMMK-------QPHVF---KVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE