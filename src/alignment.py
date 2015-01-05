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
        retrack.append((current, trace[current][0]))
        del trace[current][0]
        if len(trace[current]) == 1:
            del choice_track[-1]
            should_retrack = True
    else:
        current = (l, c)
        path = ["", ""]
    while current != (0, 0):
        all_previous = trace[current]
        previous = all_previous[0]
        direction = all_previous[0][1]
        if len(all_previous) > 1:
            op += 1
            choice_track.append([path[:], current])
        if direction == "diag":
            path[0] += s1[current[1] - 1]
            path[1] += s2[current[0] - 1]
        elif direction == "down":
            path[0] += "-"
            path[1] += s2[current[0] - 1]
        elif direction == "right":
            path[0] += s1[current[1] - 1]
            path[1] += "-"
        else:
            #free taxi ride
            pass
        current = previous[0]
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


def sequence_alignment(s1, s2):
    line = len(s2)
    column = len(s1)
    range_c = range(column + 1)
    range_l = range(line + 1)
    #trace example: {(1,2):[((0,2),"down"),((1,1),"right")]}
    #this means that the max value in node 1,2 can come from 0,2 if going down, n' from 1,1 going right
    trace = {}
    s = [[0 for _ in range_c] for _ in range_l]
    #first column and row are all 0 for now
    #column 0
    for i in range_l[1:]:
        trace.setdefault((i,0), []).append(((i - 1, 0), "down"))
    #line 0
    for i in range_c[1:]:
        trace.setdefault((0,i), []).append(((0, i - 1), "right"))

    # print "start matrix: "
    # print s
    for i in range(len(s))[1:]:
        for j in range(len(s[0]))[1:]:
            match = s2[i - 1] == s1[j - 1]
            if match:
                diagw = 1
            else:
                diagw = 0
            downward = s[i - 1][j] + 0
            rightward = s[i][j - 1] + 0
            diagonal = s[i - 1][j - 1] + diagw
            max_p = max(downward, rightward, diagonal)
            s[i][j] = max_p
            if max_p == diagonal:
                trace.setdefault((i,j), []).append(((i - 1, j - 1), "diag"))
            if max_p == rightward:
                trace.setdefault((i,j), []).append(((i,j - 1), "right"))
            if max_p == downward:
                trace.setdefault((i,j), []).append(((i - 1,j), "down"))

    print "final matrix: "
    for m in s:
        print m
    print trace
    return find_all_alignments(column, line, trace, s1, s2)


def output_alignment(s1, s2):
    alignments = sequence_alignment(s1, s2)
    print "Possible:"
    count = 0
    for a in alignments:
        count += 1
        print count, ": "
        print a[0]
        print a[1]
        print "edit distance: " + a[2]


def output_lcs(s1, s2, filew=None):
    alignments = sequence_alignment(s1, s2)
    print "Possible:"
    count = 0
    for a in alignments:
        lcs = ""
        count += 1
        if filew is None:
            print count, ": "
            for i in range(len(a[0])):
                if a[0][i] == a[1][i]:
                    lcs += a[0][i]
            print lcs
        else:
            filew.write(str(count) + ": " + "\n")
            for i in range(len(a[0])):
                if a[0][i] == a[1][i]:
                    lcs += a[0][i]
            filew.write(lcs + "\n")


def topological_order_nodes(graph, backtrace, start_nodes, source):
    aux_graph = deepcopy(graph)
    aux_bt = deepcopy(backtrace)
    # from_source = [source]
    ordered = []
    candidates = start_nodes[:]
    # candidates.remove(source)
    #todo ignore nodes that aren't a path from source
    while len(candidates) > 0:
        a = choice(candidates)
        # if a in from_source:
        ordered.append(a)
        candidates.remove(a)
        try:
            for b in graph[a].keys():
                del aux_graph[a][b]
                aux_bt[b].remove(a)
                if len(aux_bt[b]) == 0:
                    candidates.append(b)
        except KeyError:
            pass
    return ordered


def longest_path(source, sink, graph):
    # node:[(0, "source")]
    trace = {}
    #all nodes with indegree > 0
    nodes_values = {n:0 for n in set(sum([v.keys() for v in graph.values()], []))}
    start_nodes = [x for x in graph.keys() if x not in nodes_values]
    for s in start_nodes:
        if s == source:
            #maximizes value from source, so other paths wont be as good
            nodes_values[s] = 100
        else:
            nodes_values[s] = 0
    backtrace = {}
    for k, value in graph.iteritems():
        for v in value:
            backtrace.setdefault(v, []).append(k)
    ordered = topological_order_nodes(graph, backtrace, start_nodes, source)
    for n in ordered:
        if n in start_nodes:
            #no need to compute value, start nodes value is 0
            continue
        maxv = [(0, None)]
        for b in backtrace[n]:
            value = nodes_values[b] + graph[b][n]
            if value > maxv[0][0]:
                maxv = [(value, b)]
            elif value == maxv[0][0]:
                maxv.append((value, b))
        for m in maxv:
            trace.setdefault(n, []).append(m[1])
        if n == source:
            #maximizes value from source, so other paths wont be as good
            nodes_values[n] = maxv[0][0] + 100
        else:
            nodes_values[n] = maxv[0][0]
    print trace
    path = [sink]
    current = sink
    while current != source:
        all_previous = trace[current]
        prev = choice(all_previous)
        path.append(prev)
        current = prev
    return nodes_values[sink] - 100, path[::-1]


def read_longest_path():
    graph = {}
    lines = [re.split('->|:', l) for l in get_file("/data/lp_in.txt").read().splitlines()]
    for l in lines[2:]:
        graph.setdefault(l[0], {})[l[1]] = int(l[2])
    lp = longest_path(lines[0][0], lines[1][0], graph)
    print lp[0]
    print "->".join(lp[1])


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
    #this means that the max value in node 1,2 can come from 0,2 if going down, n' from 1,1 going right
    trace = {}
    if local:
        s = [[0] for _ in range_l]
        s[0] += ([0 for _ in range_c[1:]])
    elif globalt:
        s = [[0]]
        s += [[-p + (l * (-ep))] for l in range_l[:-1]]
        s[0] += ([-p + (c * (-ep)) for c in range_c[:-1]])
    else:
        #first column is 0 for fitting alignment
        s = [[0] for _ in range_l]
        s[0] += ([c * (-p) for c in range_c[1:]])

    if globalt:
        #column 0
        for i in range_l[1:]:
            trace.setdefault((i,0), []).append(((i - 1, 0), "down"))
        #line 0
        for i in range_c[1:]:
            trace.setdefault((0,i), []).append(((0, i - 1), "right"))
    elif fitting or overlap:
        #column 0
        for i in range_l[1:]:
            trace.setdefault((i,0), []).append(((i - 1, 0), "free"))
        #line 0
        for i in range_c[1:]:
            trace.setdefault((0,i), []).append(((0, i - 1), "right"))
    else:
        #column 0
        for i in range_l[1:]:
            trace.setdefault((i,0), []).append(((i - 1, 0), "free"))
            #line 0
        for i in range_c[1:]:
            trace.setdefault((0,i), []).append(((0, i - 1), "free"))

    print "start matrix: "
    for m in s:
        print m

    #for keeping the greatest value (sink) to local alignment
    maxv = -1000
    end_line = 0
    end_column = 0

    for i in range(len(s))[1:]:
        for j in range(len(s[0]))[1:]:
            sub_score = subm[s1[j - 1]][s2[i - 1]]
            diagonal = s[i - 1][j - 1] + sub_score

            down_p = p
            right_p = p
            last_if_right = trace[(i,j - 1)][0]
            if last_if_right[1] == "right":
                right_p = ep
            last_if_down = trace[(i - 1,j)][0]
            if last_if_down[1] == "down":
                down_p = ep

            downward = s[i - 1][j] - down_p
            rightward = s[i][j - 1] - right_p

            if local:
                max_p = max(downward, rightward, diagonal, 0)
            else:
                max_p = max(downward, rightward, diagonal)
            s[i].append(max_p)
            if max_p >= maxv:
                maxv = max_p
                end_line = i
                end_column = j
            if max_p == diagonal:
                trace.setdefault((i,j), []).append(((i - 1, j - 1), "diag"))
            if max_p == downward:
                trace.setdefault((i,j), []).append(((i - 1,j), "down"))
            if max_p == rightward:
                trace.setdefault((i,j), []).append(((i,j - 1), "right"))
            if local:
                if max_p == 0:
                    trace.setdefault((i,j), []).append(((0,0), "free"))

    print "final matrix: "
    for m in s:
        print m
    if globalt:
        end_line = line
        end_column = column
    elif fitting:
        max_c = -1000
        for i in range(len(s)):
            val = s[i][column]
            if val >= max_c:
                max_c = val
                end_line = i
        end_column = column
    elif overlap:
        end_line = line
        max_l = -1000
        for i in range(len(s[line])):
            val = s[line][i]
            if val >= max_l:
                max_l = val
                end_column = i
    print "FINAL SCORE: " + str(s[end_line][end_column]) + " (" + str(end_line) + "," + str(end_column) + ")"
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


# t2 = "PRTEINS"
# t1 = "PRTWPSEIN"
t2 = "YHFDVPDCWAHRYWVENPQAIAQMEQICFNWFPSMMMKQPHVFKVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE"
t1 = "YHEDVAHEDAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPIISATCARMRVRTVWE"
output_alignment_scored(t1, t2, 11, 1, "global", "blossum62")
# construct_alignment_graph(t1, t2, -5, read_subm())



#179
# YHEDV----AH----ED--AIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPI----ISATCARMRVRTVWE
# YHFDVPDCWAHRYWVENPQAIAQM-------EQICFNWFPSMMMK-------QPHVF---KVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE































def retrace_path(l, c, trace, options):
    current = (l, c)
    path = [current]
    while current != (0, 0):
        all_previous = trace[current]
        previous = all_previous[0]
        if len(all_previous) > 1:
            options += 1
            del trace[current][0]
        path.append(previous)
        current = previous
    options -= 1
    return path[::-1], options


def find_paths(column, line, trace):
    paths = []
    options = 1
    while options > 0:
        result = retrace_path(line, column, trace, options)
        paths.append(result[0])
        options = result[1]
    return paths


def alignment_graph(down, right, diag):
    line = len(down)
    column = len(right) - 1
    range_c = range(column + 1)
    range_l = range(line + 1)
    trace = {}
    s = [[0 for _ in range_c] for _ in range_l]
    #column 0
    for i in range_l[1:]:
        s[i][0] += s[i - 1][0] + down[i - 1][0]
        trace.setdefault((i,0), []).append((i - 1, 0))
    #line 0
    for i in range_c[1:]:
        s[0][i] += s[0][i - 1] + right[0][i - 1]
        trace.setdefault((0,i), []).append((0, i - 1))
    print "start matrix: "
    print s
    for i in range(line + 1)[1:]:
        for j in range(column + 1)[1:]:
            downward = s[i - 1][j] + down[i - 1][j]
            rightward = s[i][j - 1] + right[i][j - 1]
            diagonal = s[i - 1][j - 1] + diag[i - 1][j - 1]
            max_p = max(downward, rightward, diagonal)
            s[i][j] = max_p
            if max_p == diagonal:
                trace.setdefault((i,j), []).append((i - 1, j - 1))
            if max_p == downward:
                trace.setdefault((i,j), []).append((i - 1,j))
            if max_p == rightward:
                trace.setdefault((i,j), []).append((i,j - 1))
    print "final matrix: "
    print s
    print trace
    return find_paths(column, line, trace)


def read_matrixes(filename):
    data = get_file(filename).read().splitlines()
    first_split = data.index("-")
    down = [[int(x) for x in d.split(" ")] for d in data[:first_split]]
    print "down: "
    print down
    last_split = -(data[::-1].index("-"))
    right = [[int(x) for x in d.split(" ")] for d in data[first_split + 1:last_split - 1]]
    print "right: "
    print right
    diag = [[int(x) for x in d.split(" ")] for d in data[last_split:]]
    print "diag: "
    print diag
    return alignment_graph(down, right, diag)


def construct_alignment_graph(s1, s2, p, subm):
    graph = {}
    for x in s1:
        for y in s2:
            # graph.setdefault(x, {})[y] = subm[]
            return

def construct_alignment(l, c, trace, s1, s2, op):
    current = (l, c)
    path = ["", ""]
    deletedpath = False
    while current != (0, 0):
        all_previous = trace[current]
        previous = all_previous[0]
        direction = all_previous[0][1]
        if len(all_previous) > 1:
            op += 1
            if not deletedpath:
                del trace[current][0]
                deletedpath = True
        if direction == "diag":
            path[0] += s1[current[1] - 1]
            path[1] += s2[current[0] - 1]
        elif direction == "down":
            path[0] += "-"
            path[1] += s2[current[0] - 1]
        elif direction == "right":
            path[0] += s1[current[1] - 1]
            path[1] += "-"
        else:
            #free taxi ride
            pass
        current = previous[0]

    distance = hamming_d(path[0], path[1])
    return [path[0][::-1], path[1][::-1], distance], op