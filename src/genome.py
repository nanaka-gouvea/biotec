from random import choice

__author__ = 'natalia'
from translation import get_file
from translation import get_file_w
from collections import defaultdict
from copy import deepcopy


def composition(k, text):
    comp = []
    for i in range(len(text) - k + 1):
        comp.append(text[i:i + k])
    return sorted(comp)


def string_spelled(peaces):
    k = len(peaces[0])
    genome = peaces[0][:k - 1]
    for p in peaces:
        genome += p[-1]
    return genome


def overlap_graph(peaces):
    dir_g = dict.fromkeys(peaces, [])
    for suffix in peaces:
        prefixes = []
        for prefix in peaces:
            if suffix[1:] == prefix[:-1]:
                prefixes.append(prefix)
        dir_g[suffix] = prefixes
    return dir_g


def de_bruijn(k, text):
    peaces = [text[i:i + k] for i in range(len(text) - k + 1)]
    nodes = [peaces[0][:k - 1]]
    for p in peaces:
        nodes.append(p[1:])
    graph = {n: [] for n in set(nodes)}
    for i in range(len(nodes) - 1):
        graph[nodes[i]].append(nodes[i + 1])
    return graph


def de_bruijn_patterns(peaces):
    peaces = sorted(peaces)
    k = len(peaces[0])
    unordered_edges = []
    edge_counter = defaultdict(int)
    for p in peaces:
        edge = (p[:k - 1], p[1:])
        unordered_edges.append(edge)
        edge_counter[edge] += 1
    graph = {}
    for e in unordered_edges:
        if edge_counter[e] > 0:
            graph.setdefault(e[0], []).append(e[1])
            edge_counter[e] -= 1
        aux = unordered_edges[:]
        aux.remove(e)
        for next_e in aux:
            if e[1] == next_e[0] and edge_counter[next_e] > 0:
                graph.setdefault(e[1], []).append(next_e[1])
                edge_counter[next_e] -= 1
    return graph


def output_graph(graph, filew=None):
    for peace in sorted(graph.keys()):
        endings = graph[peace]
        if len(endings) > 0:
            if filew is None:
                print peace, "->", ','.join(endings)
            else:
                filew.write(''.join([peace, " -> ", ','.join(endings), "\n"]))


def remove_step(new_step, walking_path):
    for k, v in walking_path.iteritems():
        if new_step in v:
            walking_path[k].remove(new_step)


def eulerian_cycle(graph):
    last_cycle = []
    cycle = []
    nodes = graph.keys()
    walking_path = deepcopy(graph)
    l_edges = len(sum(graph.values(), []))
    while len(cycle) < l_edges:
        cycle = []
        if len(last_cycle) == 0:
            cycle.append(choice(nodes))
        else:
            # last_path_cycle = last_cycle[:-1]
            unexplored = 0
            for i in range(len(last_cycle)):
                if len(walking_path[last_cycle[i]]) > 0:
                    unexplored = i
                    break
            # transverse
            cycle = last_cycle[unexplored:] + last_cycle[1:unexplored + 1]
        step = cycle[-1]
        while len(walking_path[step]) > 0:
            new_step = choice(walking_path[step])
            cycle.append(new_step)
            # o novo no nao podera ser visitado novamente nesse cycle
            walking_path[step].remove(new_step)
            step = new_step
        last_cycle = cycle[:]
    return cycle


def eulerian_path(graph):
    # balance graph
    total_ins = sum(graph.values(), [])
    miss_out = None
    miss_in = None
    nodes = set(sum(graph.values(), graph.keys()))
    for n in nodes:
        try:
            v = graph[n]
            out = len(v)
        except KeyError:
            out = 0
        ins = total_ins.count(n)
        if out < ins:
            miss_out = n
        elif out > ins:
            miss_in = n
    if miss_in is not None and miss_out is not None:
        try:
            graph[miss_out].append(miss_in)
        except KeyError:
            graph[miss_out] = [miss_in]

        path = eulerian_cycle(graph)
        #split at added edge
        spliti = 0
        for i in range(len(path) - 1):
            if path[i] == miss_out and path[i + 1] == miss_in:
                spliti = i + 1
        path = path[spliti:] + path[1:spliti]
        return path
    else:
        return eulerian_cycle(graph)[:-1]


# sample = get_file("/data/euler_in.txt").read().splitlines()
# smap = {}
# for line in sample:
#     parts = line.split(" -> ")
#     smap[parts[0]] = parts[1].split(",")
# print smap
# print "->".join(eulerian_path(smap))
# fileo = get_file_w("/data/euler_out.txt")
# fileo.write("->".join(eulerian_path(smap)))

sample = "CTTA ACCA TACC GGCT GCTT TTAC".split(" ")
dbr = de_bruijn_patterns(sample)
print dbr
eu = eulerian_path(dbr)
print eu
print eu[0] + "".join(eu[1:])









# TODO acertar isso ae
def kuniversal_binary(k):
    i = 0
    ku = "000"
    while True:
        b = bin(i)[2:]
        lb = len(b)
        if lb == k and b not in ku:
            ku += b
        elif lb > k:
            break
        i += 1
    return ku

