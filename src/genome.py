__author__ = 'natalia'
from translation import get_file
from translation import get_file_w


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
    k = len(peaces[0])
    dir_g = dict.fromkeys(peaces,[])
    for ps in peaces:
        for pp in peaces:
            if ps[1:] == pp[:-1]:
                dir_g[ps].append(pp)
    return dir_g


sample = "ATGCG GCATG CATGC AGGCA GGCAT".split(" ")
print overlap_graph(sample)

# teste = {"TESTEA":[], "TESTEB":[]}
# teste["TESTEA"].append("VALED")
# print teste