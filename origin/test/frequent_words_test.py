from unittest.case import TestCase
from origin.main import frequent_words

__author__ = 'natalia'
import unittest

class FrequentWordsTest(TestCase):
    def test_clump_3_3_1(self):
        clumps = frequent_words.find_clump_patterns("BAG", 3, 3, 1)
        assert len(clumps) == 1
        assert "BAG" in clumps