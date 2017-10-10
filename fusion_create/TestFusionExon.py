# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 14:11:13 2017

@author: aelamb
"""

from FusionExon import FusionExon as FEx
import unittest
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import random


class TestFusionEvent(unittest.TestCase):
    
    def setUp(self):
        self.exon1 = SeqFeature(FeatureLocation(ExactPosition(100), 
                                                ExactPosition(500), 
                                                strand = -1),
                          type = 'exon', id = 'exon_1')
        self.exon2 = SeqFeature(FeatureLocation(ExactPosition(800), 
                                                ExactPosition(1000), 
                                                strand = -1),
                          type = 'exon', id = 'exon_2')
        self.exon3 = SeqFeature(FeatureLocation(ExactPosition(1500), 
                                                ExactPosition(2000), 
                                                strand = -1),
                          type = 'exon', id = 'exon_3')
        self.FEx1 = FEx(self.exon1)
                          
    def test_getters(self):
        self.assertEqual(self.FEx1.get_start(), 100)
        self.assertEqual(self.FEx1.get_end(), 500)
    
    def test_setters(self):
        self.FEx1.set_start(150)
        self.FEx1.set_end(200)
        self.assertEqual(self.FEx1.get_start(), 150)
        self.assertEqual(self.FEx1.get_end(), 200)
        self.setUp()
    
    def test_is_breakable(self):
        self.assertTrue(self.FEx1.is_breakable(300, 100))
        self.assertFalse(self.FEx1.is_breakable(300, 101))
    
    def test_break_exon_at_junction(self):
        random.seed(1)
        FEx = self.FEx1
        self.assertEqual(FEx.break_exon_at_junction(1, 1, "forward"), 154)
        self.assertEqual(FEx.get_start(), 100)
        self.assertEqual(FEx.get_end(), 154)
        self.setUp()
        FEx = self.FEx1
        self.assertEqual(FEx.break_exon_at_junction(1, 1, "reverse"), 439)
        self.assertEqual(FEx.get_start(), 439)
        self.assertEqual(FEx.get_end(), 500)
        self.setUp()

if __name__ == '__main__':
    unittest.main()