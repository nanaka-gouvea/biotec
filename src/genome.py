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


def string_spelled(pieces):
    k = len(pieces[0])
    genome = pieces[0][:k - 1]
    for p in pieces:
        genome += p[-1]
    return genome


def overlap_graph(pieces):
    dir_g = dict.fromkeys(pieces, [])
    for suffix in pieces:
        prefixes = []
        for prefix in pieces:
            if suffix[1:] == prefix[:-1]:
                prefixes.append(prefix)
        dir_g[suffix] = prefixes
    return dir_g


def de_bruijn(k, text):
    pieces = [text[i:i + k] for i in range(len(text) - k + 1)]
    nodes = [pieces[0][:k - 1]]
    for p in pieces:
        nodes.append(p[1:])
    graph = {n: [] for n in set(nodes)}
    for i in range(len(nodes) - 1):
        graph[nodes[i]].append(nodes[i + 1])
    return graph


def de_bruijn_patterns(pieces):
    pieces = sorted(pieces)
    k = len(pieces[0])
    unordered_edges = []
    edge_counter = defaultdict(int)
    for p in pieces:
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
    for piece in sorted(graph.keys()):
        endings = graph[piece]
        if len(endings) > 0:
            if filew is None:
                print piece, "->", ','.join(endings)
            else:
                filew.write(''.join([piece, " -> ", ','.join(endings), "\n"]))


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
    while len(cycle) <= l_edges:
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


def reconstruction(reads):
    k = len(reads[0])
    dbr = de_bruijn_patterns(reads)
    eu = eulerian_path(dbr)
    return eu[0] + "".join([x[k - 2:] for x in eu[1:]])


def generate_binary_strings(lenght):
    max_b = sum([2 ** x for x in range(lenght)])
    return [bin(x)[2:].zfill(lenght) for x in range(max_b + 1)]


def universal_circular_string(k):
    pieces = generate_binary_strings(k)
    print pieces
    dbr = de_bruijn_patterns(pieces)
    print dbr
    eu = eulerian_cycle(dbr)
    print eu
    return "".join([x[0] for x in eu[:-1]])


def universal_string(k):
    pieces = generate_binary_strings(k)
    print pieces
    dbr = de_bruijn_patterns(pieces)
    print dbr
    eu = eulerian_path(dbr)
    print eu
    return "".join([x[0] for x in eu])


def kd_mer_composition(k, d, string):
    comp = []
    for i in range(len(string) - (k + d) - (k - 1)):
        kmer1 = string[i:i + k]
        kmer2 = string[i + k + d:i + k + d + k]
        comp.append((kmer1, kmer2))
    return sorted(comp)

for c in kd_mer_composition(3,2,"TAATGCCATGGGATGTT"):
    print "(" + c[0] + "|" + c[1] + ")"




# sample = get_file("/data/euler_in.txt").read().splitlines()
