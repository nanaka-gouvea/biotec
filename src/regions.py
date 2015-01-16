__author__ = 'natalia'
from translation import get_file_w
from translation import get_file
from itertools import cycle


def greedy_sorting(p, out=None):
    aprox_distance = 0
    for k in range(len(p) + 1)[1:]:
        i = k - 1
        if p[i] != k:
            index_of_should_be = [abs(n) for n in p].index(k)
            reverse_block = [-1 * n for n in p[i:index_of_should_be + 1]]
            reverse_block.reverse()
            p = p[:i] + reverse_block + p[index_of_should_be + 1:]
            aprox_distance += 1
            if out is None:
                print "(" + " ".join("+" + str(n) if n > 0 else str(n) for n in p) + ")"
            else:
                out.write("(" + " ".join("+" + str(n) if n > 0 else str(n) for n in p) + ")\n")
        if p[i] == -k:
            #k-sorting reversal
            p[i] *= -1
            aprox_distance += 1
            if out is None:
                print "(" + " ".join("+" + str(n) if n > 0 else str(n) for n in p) + ")"
            else:
                out.write("(" + " ".join("+" + str(n) if n > 0 else str(n) for n in p) + ")\n")
    print aprox_distance
    return aprox_distance


def count_breaks(p):
    breaks = 0
    p.append(len(p) + 1)
    p = [0] + p
    for i in range(len(p) - 1):
        a = p[i]
        b = p[i + 1]
        if b != a + 1:
            breaks += 1
    return breaks


def count_permutation(p):
    perm = 0
    for i in range(len(p) - 1):
        a = p[i]
        b = p[i+1]
        if b == a + 1:
            perm += 1


def chromosome_to_cycle(chromosome):
    nodes = [None for _ in range(len(chromosome)*2)]
    for j in range(len(chromosome)):
        i = chromosome[j]
        if i > 0:
            nodes[2*j] = 2*i - 1
            nodes[2*j + 1] = 2*i
        else:
            nodes[2*j] = -2*i
            nodes[2*j + 1] = -2*i -1
    return nodes


def cycle_to_chromossome(nodes):
    chromossome = [0 for n in range(len(nodes)/2)]
    for j in range(len(nodes)/2):
        if nodes[2*j] < nodes[2*j + 1]:
            chromossome[j] = nodes[2*j + 1] / 2
        else:
            chromossome[j] = -nodes[2*j] / 2
    return chromossome


def colored_edges(p):
    edges = []
    for c in p:
        nodes = chromosome_to_cycle(c)
        nodes.append(nodes[0])
        for j in range(len(c)):
            edges.append((nodes[2*j+1], nodes[2*j+2]))
    return edges


#XXX what?
def graph_to_genome(graph):
    p = []
    for n in graph:
        c = cycle_to_chromossome(n)
        p.append(c)
    return p


# greedy_sorting([+20, +7, +10, +9, +11, +13, +18, -8, -6, -14, +2, -4, -16, +15, +1, +17, +12, -5, +3, -19])
# print count_breaks([-16, -20, +11, +12, -14, -13, -15, -6, -8, -19, -18, -17, -10, +4, -5, -2, +7, -3, +1, -9])











