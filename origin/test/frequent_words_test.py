from unittest.case import TestCase

from origin import frequent_words


__author__ = 'natalia'


class FrequentWordsTest(TestCase):
    def test_clump_3_3_1(self):
        clumps = frequent_words.find_clump_patterns("BAG", 3, 3, 1)
        assert len(clumps) == 1
        assert "BAG" in clumps

    def test_frequent_small(self):
        f = frequent_words.most_frequent_words(open("../data/ecoli_region.txt").read(), 9)
        print f