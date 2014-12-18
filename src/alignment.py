from sys import maxint
from translation import get_file
from translation import get_file_w
from genome import non_branching_paths

__author__ = 'natalia'


def substitution_matrix():
    subst_matrix = {}
    subst_matrix["A"] = {"C":-5, "T":-5, "G":-7, "A":2}
    subst_matrix["C"] = {"C":2, "T":-7, "G":-5, "A":-5}
    subst_matrix["G"] = {"C":-5, "T":-7, "G":2, "A":-7}
    subst_matrix["T"] = {"C":-7, "T":2, "G":-7, "A":-5}
    return subst_matrix


def construct_alignment(l, c, trace, s1, s2, options):
    current = (l, c)
    path = ["", ""]
    while current != (0, 0):
        all_previous = trace[current]
        previous = all_previous[0]
        direction = all_previous[0][1]
        if len(all_previous) > 1:
            options += 1
            del trace[current][0]
        if direction == "diag":
            path[0] += s1[current[1] - 1]
            path[1] += s2[current[0] - 1]
        elif direction == "down":
            path[0] += "-"
            path[1] += s2[current[0] - 1]
        else:
            path[0] += s1[current[1] - 1]
            path[1] += "-"
        current = previous[0]
    options -= 1
    return path[::-1], options


def find_all_alignments(column, line, trace, s1, s2):
    paths = []
    options = 1
    while options > 0:
        result = construct_alignment(line, column, trace, s1, s2, options)
        paths.append(result[0])
        options = result[1]
    return paths


def sequence_alignment(s1, s2):
    line = len(s2)
    column = len(s1)
    range_c = range(column + 1)
    range_l = range(line + 1)
    trace = {}
    s = [[0 for _ in range_c] for _ in range_l]
    #first column and row are all 0 for now
    #column 0
    for i in range_l[1:]:
        trace.setdefault((i,0), []).append(((i - 1, 0), "down"))
    #line 0
    for i in range_c[1:]:
        trace.setdefault((0,i), []).append(((0, i - 1), "right"))

    print "start matrix: "
    print s
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
            if max_p == downward:
                trace.setdefault((i,j), []).append(((i - 1,j), "down"))
            if max_p == rightward:
                trace.setdefault((i,j), []).append(((i,j - 1), "right"))
    print "final matrix: "
    print s
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

# t1 = "AACCTTGG"
# t2 = "ACACTGTGA"
# output_alignment(t1, t2)


def output_lcs(backtrack, v, i, j):
    if i == 0 or j == 0:
        return
    if backtrack[i][j] == "down":
        output_lcs(backtrack, v, i - 1, j)
    elif backtrack[i,j] == "right":
        output_lcs(backtrack, v, i, j - 1)
    else:
        output_lcs(backtrack, v, i - 1, j - 1)
        print v[i]














































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