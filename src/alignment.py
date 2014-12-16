from sys import maxint
from itertools import cycle

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



# def dp_change(money, coins):
#     min_num_coins = {0:0}
#     for m in range(1, money):
#         min_num_coins[m] = maxint
#         for c in coins:
#             if m >= c:
#                 if dpc_change()


# print min_num_coins_change(22, [5,4,1])
# print min_num_coins_change(16053, [23,22,21,14,6,5,3,1])
