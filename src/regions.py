__author__ = 'natalia'
from translation import get_file_w
from translation import get_file


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


# greedy_sorting([-3, +4, +1, +5, -2], get_file_w("/data/greedy_sort_out.txt"))
# get_file_w("/data/greedy_sort_out.txt").write("teste")

# print count_breaks([+3, +4, +5, -12, -8, -7, -6, +1, +2, +10, +9, -11, +13, +14])

permut = []
for n in get_file("/data/number_of_breaks.txt").readline().split(" "):
    permut.append(int(n))
print count_breaks(permut)