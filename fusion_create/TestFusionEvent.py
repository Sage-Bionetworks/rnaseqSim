# -*- coding: utf-8 -*-

from FusionEvent import FusionEvent
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import unittest
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
        self.exon4 = SeqFeature(FeatureLocation(ExactPosition(10000), 
                                                ExactPosition(15000), 
                                                strand = -1),
                          type = 'exon', id = 'exon_4')
        self.exon5 = SeqFeature(FeatureLocation(ExactPosition(17000), 
                                                ExactPosition(20000), 
                                                strand = -1),
                          type = 'exon', id = 'exon_5')
        self.exon6 = SeqFeature(FeatureLocation(ExactPosition(22000), 
                                                ExactPosition(25000), 
                                                strand = -1),
                          type = 'exon', id = 'exon_6')
        self.exons1 = [self.exon1, self.exon2, self.exon3]
        self.exons2 = [self.exon4, self.exon5, self.exon6]
        self.FE1 = FusionEvent(self.exons1, self.exons2, "-", "-")
        self.FE2 = FusionEvent(self.exons1, self.exons2, "+", "+")
                          

    def test_getters(self):
        self.assertEqual(self.FE1.get_acceptor_exons(), self.exons2)
        self.assertEqual(self.FE1.get_donor_exons(), self.exons1)
        self.assertEqual(self.FE1.get_acceptor_strand(), "-")
        self.assertEqual(self.FE1.get_donor_strand(), "-")
        self.assertEqual(self.FE1.get_acceptor_junction(), None)
        self.assertEqual(self.FE1.get_donor_junction(), None)
    
    def test_create_breakage(self):
        # both strands are "-", donor gets cut off on left, acceptor on right
        random.seed(1)
        res1 = self.FE1.create_breakage("donor")
        self.assertEqual(res1[0], 799)
        self.assertEqual(res1[1], [self.exon2, self.exon3])
        random.seed(1)
        res2 = self.FE1.create_breakage("acceptor")
        self.assertEqual(res2[0], 14999)
        self.assertEqual(res2[1], [self.exon4])
        # both strands are "+", donor gets cut off on right, acceptor on left
        random.seed(1)
        res3 = self.FE2.create_breakage("donor")
        self.assertEqual(res3[0], 499)
        self.assertEqual(res3[1], [self.exon1])
        random.seed(4)
        res4 = self.FE2.create_breakage("acceptor")
        self.assertEqual(res4[0], 16999)
        self.assertEqual(res4[1], [self.exon5, self.exon6])
        
        # mid exon breakages
        random.seed(2)
        res5 = self.FE1.create_breakage("donor", True)
        self.assertEqual(res5[0], 1972)
        self.assertEqual(res5[1][0].location.start, ExactPosition(1973))
        self.assertEqual(res5[1][0].location.end, ExactPosition(2000))
        random.seed(2)
        res6 = self.FE1.create_breakage("acceptor", True)
        self.assertEqual(res6[0], 24842)
        self.assertEqual(res6[1][0].location.start, ExactPosition(10000))
        self.assertEqual(res6[1][0].location.end, ExactPosition(15000))
        self.assertEqual(res6[1][1].location.start, ExactPosition(17000))
        self.assertEqual(res6[1][1].location.end, ExactPosition(20000))                 
        self.assertEqual(res6[1][2].location.start, ExactPosition(22000))
        self.assertEqual(res6[1][2].location.end, ExactPosition(24843))                             

if __name__ == '__main__':
    unittest.main()