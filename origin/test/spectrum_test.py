from unittest.case import TestCase
from origin.spectrum import *
__author__ = 'natalia'


class Spectrum_Test(TestCase):

    def test_theoretical_spectrum(self):
        lines = open("../data/theoretical_spectrum.txt").read().splitlines()
        ts = theoretical_spectrum(lines[0])
        self.assertEqual(ts, [int(x) for x in lines[1].split(" ")])


