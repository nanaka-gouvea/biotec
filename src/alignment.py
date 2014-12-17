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


def global_alignment(s1, s2, d):
    if len(s1) != len(s2):
        return
    # am = [[] for _ in "-" + s1]
    ps = 0
    # px = 0
    # py = 0
    dseq = range(0, d * len(s1) - 1, d)
    am = [dseq]
    for _ in range(len("-" + s1)):
        am.append("a")
    # for a in am:


# sm = substitution_matrix()
# print sm["A"]["G"]


def min_num_coins_change(money, coins, prev=None):
    if money == 0:
        return 0
    if prev is None:
        prev = {}
    min_num_coins = maxint
    for c in coins:
        if money >= c:
            prevm = money - c
            try:
                num_coins = prev[prevm]
            except KeyError:
                num_coins = min_num_coins_change(prevm, coins, prev)
                prev[prevm] = num_coins
                if len(prev.keys()) > coins[0]:
                    del prev[min(prev.keys())]
            if num_coins + 1 < min_num_coins:
                min_num_coins = num_coins + 1
    return min_num_coins


def change(money, coins, prev=None):
    if money == 0:
        return []
    if prev is None:
        prev = {}
    min_coins = range(money)
    for c in coins:
        if money >= c:
            prevm = money - c
            try:
                num_coins = prev[prevm]
            except KeyError:
                num_coins = change(prevm, coins, prev)
                prev[prevm] = num_coins
                if len(prev.keys()) > coins[0]:
                    del prev[min(prev.keys())]
            if len(num_coins) + 1 < len(min_coins):
                #todo nao adicionar '0'
                min_coins = num_coins
                min_coins.append(c)
    return min_coins


def manhattam_tourist(line, column, down, right):
    """
        down and right = [[]] (0,1) = down[0][1]
    """
    range_m = range(column + 1)
    range_n = range(line + 1)
    s = [[0 for _ in range_m] for _ in range_n]
    #column 0
    for i in range_n[1:]:
        s[i][0] += s[i - 1][0] + down[i - 1][0]
    #line 0
    for i in range_m[1:]:
        s[0][i] += s[0][i - 1] + right[0][i - 1]
    print "start matrix: "
    print s
    for i in range(line + 1)[1:]:
        for j in range(column + 1)[1:]:
            downward = s[i - 1][j] + down[i - 1][j]
            rightward = s[i][j - 1] + right[i][j - 1]
            s[i][j] = max(downward, rightward)
    print "final matrix: "
    print s
    return s[line][column]


def read_manhattam_tourist(filename):
    data = get_file(filename).read().splitlines()
    n = int(data[0].split(" ")[0])
    m = int(data[0].split(" ")[1])
    split_i = data.index("-")
    down = [[int(x) for x in d.split(" ")] for d in data[1:split_i]]
    print "down: "
    print down
    right = [[int(x) for x in d.split(" ")] for d in data[split_i + 1:]]
    print "right: "
    print right
    return manhattam_tourist(n, m, down, right)


def retrace_path(l, c, d_diag, d_down, trace, options):
    path = []
    current = (l, c)
    while current != (0, 0):
        directions = trace[current]
        direction = directions[0]
        # todo tratar mais de um alinhamento otimo
        if len(directions) > 1:
            options += 1
            del trace[current][0]
        path.append(direction)
        if direction == d_diag:
            current = (current[0] - 1, current[1] - 1)
        elif direction == d_down:
            current = (current[0] - 1, current[1])
        else:
            current = (current[0], current[1] - 1)
    options -= 1
    return path[::-1], options


def best_path(down, right, diag):
    d_diag = "s"
    d_down = "1"
    d_right = "2"
    line = len(down)
    column = len(right) - 1
    range_c = range(column + 1)
    range_l = range(line + 1)
    s = [[0 for _ in range_c] for _ in range_l]
    #column 0
    for i in range_l[1:]:
        s[i][0] += s[i - 1][0] + down[i - 1][0]
    #line 0
    for i in range_c[1:]:
        s[0][i] += s[0][i - 1] + right[0][i - 1]
    print "start matrix: "
    print s
    trace = {}
    for i in range(line + 1)[1:]:
        for j in range(column + 1)[1:]:
            downward = s[i - 1][j] + down[i - 1][j]
            rightward = s[i][j - 1] + right[i][j - 1]
            diagonal = s[i - 1][j - 1] + diag[i - 1][j - 1]
            max_p = max(downward, rightward, diagonal)
            s[i][j] = max_p
            from_ = []
            if max_p == diagonal:
                #s for substitution
                from_.append(d_diag)
            if max_p == downward:
                #1 for gap in first
                from_.append(d_down)
            if max_p == rightward:
                #2 for gap in second
                from_.append(d_right)
            trace[(i,j)] = from_
    print "final matrix: "
    print s
    print trace
    paths = []
    options = 1
    while options > 0:
        result = retrace_path(line, column, d_diag, d_down, trace, options)
        paths.append(result[0])
        options = result[1]
    return paths


def best_path_graph(down, right, diag):
    line = len(down)
    column = len(right) - 1
    range_c = range(column + 1)
    range_l = range(line + 1)
    s = [[0 for _ in range_c] for _ in range_l]
    #column 0
    for i in range_l[1:]:
        s[i][0] += s[i - 1][0] + down[i - 1][0]
    #line 0
    for i in range_c[1:]:
        s[0][i] += s[0][i - 1] + right[0][i - 1]
    print "start matrix: "
    print s
    trace = {}
    for i in range(line + 1)[1:]:
        for j in range(column + 1)[1:]:
            downward = s[i - 1][j] + down[i - 1][j]
            rightward = s[i][j - 1] + right[i][j - 1]
            diagonal = s[i - 1][j - 1] + diag[i - 1][j - 1]
            max_p = max(downward, rightward, diagonal)
            s[i][j] = max_p
            if max_p == diagonal:
                #s for substitution
                trace.setdefault((i - 1, j - 1), []).append((i,j))
            if max_p == downward:
                #1 for gap in first
                trace.setdefault((i - 1,j), []).append((i,j))
            if max_p == rightward:
                #2 for gap in second
                trace.setdefault((i,j - 1), []).append((i,j))
    print "final matrix: "
    print s
    print trace
    paths = []
    nbs = non_branching_paths(trace)
    for nb in nbs:
        if nb[0] == (0,0) or nb[-1] == (line, column):
            paths.append(nb)
    return paths


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
    return best_path(down, right, diag)
    # return best_path_graph(down, right, diag)

# print read_manhattam_tourist("/data/manhattam_in.txt")
print read_matrixes("/data/manhattam_diag_in.txt")

# print min_num_coins_change(22, [5,4,1])
# print change(22, [5,4,1])
# print min_num_coins_change(16053, [23,22,21,14,6,5,3,1])