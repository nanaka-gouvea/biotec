# coding=utf-8
import numpy as np
class SubstitutionMatrix:
    'Representa uma matriz de substituição de resíduos
# TODO criar construtor q recebe a lista de pares com score e monta a matriz para cada [(x,y), score], montar [x[y]] = score & [y[x]] =score'

# esses vao para o actor
# nosub = 3
# gs = -5
# ge = 1

# subm = np.zeros(5)
# subm = [[0 for x in xrange(5)] for x in xrange(5)]
# subm['A']['C'] = -5
# subm = {'A':{'C':-5}, 'C':{'A'}}

subm = {}
subm['A'] = {'C':-5, 'T':-7, 'G':-8}
subm['C'] = {'A':-5, 'T':-3, 'G':-4}
subm['T'] = {'A':-7, 'C':-3, 'G':-6}
subm['G'] = {'A':-8, 'C':-4, 'T':-7}


def subScore(orig, sub)