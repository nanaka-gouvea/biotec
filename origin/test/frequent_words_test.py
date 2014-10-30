from unittest.case import TestCase

from origin import frequent_words


__author__ = 'natalia'


class FrequentWordsTest(TestCase):

    def test_clump_edge(self):
        clumps = frequent_words.clump_finding("TAG", 3, 3, 1)
        self.assertEqual(clumps, set(['TAG']))

        clumps = frequent_words.clump_finding("AGAGA", 3, 4, 2)
        assert len(clumps) == 0

        clumps = frequent_words.clump_finding("AGAG", 1, 4, 2)
        self.assertEqual(clumps, set(['A', 'G']))

    def test_clump_overlap(self):
        clumps = frequent_words.clump_finding("AAAA", 3, 4, 2)
        self.assertEqual(clumps, set(['AAA']))

        clumps = frequent_words.clump_finding("GCTGGCTGGCTGGCTGG", 9, 17, 3)
        self.assertEqual(clumps, set(['GCTGGCTGG']))

        clumps = frequent_words.clump_finding("AGAGAGAGAGA", 3, 9, 4)
        self.assertEqual(clumps, set(['GAG', 'AGA']))

    def test_clump_many_matches(self):
        clumps = frequent_words.clump_finding("AAAATTATATAAAAAAAAAAATTTATTTTTTTTTTATTTTTT", 3, 10, 3)
        self.assertEqual(clumps, set(['AAA', 'TTT']))

    def test_clump_medium(self):
        clumps = frequent_words.clump_finding(open("../data/test_clump.txt").read(), 12, 595, 19)
        self.assertEqual(clumps, set(['AGAGTGATTGCG', 'GTGGATAGCCTA', 'GTGATCCACCGA', 'GATAGTTGGTCT', 'ACTTCCAAACAG',
                                      'TACTCCTGAAGT', 'TTGCAAACTGAC', 'CCGCACGAAGTA', 'ATAACGATTTCC']))

    def test_clump_medium_not_over(self):
        clumps = frequent_words.clump_finding(open("../data/test_clump.txt").read(), 12, 595, 19)