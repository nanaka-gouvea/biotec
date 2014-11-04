from unittest.case import TestCase
from origin import aminoacids
__author__ = 'natalia'

class AminoacidsTest(TestCase):

    def test_clump_edge(self):
        lines = open("../data/translate.txt").read().splitlines()
        trans = aminoacids.translate_peptide(lines[0])
        self.assertEqual(trans, lines[1])