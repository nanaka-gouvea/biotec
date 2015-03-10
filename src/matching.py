__author__ = 'natalia'

def trie(patterns):
    root = "-"
    trie = {root:[]}
    for p in patterns:
        current = root
        for c in p:
            if current == c:
                current = c
            else:
                trie[current].append()

    return 0