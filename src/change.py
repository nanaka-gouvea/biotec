__author__ = 'natalia'

from sys import maxint
from translation import get_file
from translation import get_file_w


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