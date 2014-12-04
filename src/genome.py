__author__ = 'natalia'
from translation import get_file
from translation import get_file_w
from collections import defaultdict

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
    dir_g = dict.fromkeys(peaces,[])
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
    graph = {n:[] for n in set(nodes)}
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


sample = get_file("/data/deb_test.txt").read().splitlines()
# output_graph(de_bruijn_patterns(sample))
fileo = get_file_w("/data/deb_result.txt")
output_graph(de_bruijn_patterns(sample), fileo)


# expected = get_file("/data/deb_test.txt").read().splitlines()[1123:]
# result = get_file("/data/deb_result.txt").read().splitlines()
# print expected == result








#TODO acertar isso ae
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

