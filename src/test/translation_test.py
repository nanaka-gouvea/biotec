from unittest.case import TestCase
from src import translation
__author__ = 'natalia'


class TranslationTest(TestCase):

    def test_translate(self):
        lines = open("../data/translate.txt").read().splitlines()
        trans = translation.translate_peptide(lines[0])
        self.assertEqual(trans.replace(" ", ""), lines[1])

    def test_find_peptide_encoding(self):
        lines = open("../data/peptide_encoding.txt").read().splitlines()
        strings = translation.find_peptide_encoding(lines[0], lines[1])
        self.assertEqual(set(strings), set(lines[2:]))

    def test_better_find_peptide_encoding(self):
        lines = open("../data/peptide_encoding.txt").read().splitlines()
        result = translation.better_find_peptide_encoding(lines[0], lines[1])
        self.assertEqual(len(lines[2:]), len(result))
        self.assertEqual(set(lines[2:]), set(result))

