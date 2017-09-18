# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 09:33:08 2017

@author: aelamb
"""


from FusionExonList import FusionExonList
from FusionExon import FusionExon
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import unittest
import random

class TestFusionEventList(unittest.TestCase):
    
    def setUp(self):
        self.sf1 = SeqFeature(FeatureLocation(ExactPosition(100), 
                                                ExactPosition(500), 
                                                strand = -1),
                          type = 'exon', id = 'exon_1')
        self.sf2 = SeqFeature(FeatureLocation(ExactPosition(800), 
                                                ExactPosition(1000), 
                                                strand = -1),
                          type = 'exon', id = 'exon_2')
        self.sf3 = SeqFeature(FeatureLocation(ExactPosition(1500), 
                                                ExactPosition(2000), 
                                                strand = -1),
                          type = 'exon', id = 'exon_3')
        self.sf4 = SeqFeature(FeatureLocation(ExactPosition(10000), 
                                                ExactPosition(15000), 
                                                strand = -1),
                          type = 'exon', id = 'exon_4')
        self.sf5 = SeqFeature(FeatureLocation(ExactPosition(17000), 
                                                ExactPosition(20000), 
                                                strand = -1),
                          type = 'exon', id = 'exon_5')
        self.sf6 = SeqFeature(FeatureLocation(ExactPosition(22000), 
                                                ExactPosition(25000), 
                                                strand = -1),
                          type = 'exon', id = 'exon_6')
                          
        self.exon1 = FusionExon(self.sf1)
        self.exon2 = FusionExon(self.sf2)
        self.exon3 = FusionExon(self.sf3)
        self.exon4 = FusionExon(self.sf4)
        self.exon5 = FusionExon(self.sf5)
        self.exon6 = FusionExon(self.sf6)

        self.sfs1 = [self.sf1, self.sf2, self.sf3]
        self.sfs2 = [self.sf4, self.sf5, self.sf6]
        self.exons1 = [self.exon1, self.exon2, self.exon3]
        self.exons2 = [self.exon4, self.exon5, self.exon6]
        self.EXL1 = FusionExonList(self.sfs1, "+", "donor")
        self.EXL2 = FusionExonList(self.sfs1, "-", "acceptor")
        self.EXL3 = FusionExonList(self.sfs1, "+", "acceptor")
        self.EXL4 = FusionExonList(self.sfs1, "-", "donor")
    
    def test_set_direction(self):
        self.assertEqual(self.EXL1.direction, "forward")
        self.assertEqual(self.EXL2.direction, "forward")
        self.assertEqual(self.EXL3.direction, "reverse")
        self.assertEqual(self.EXL4.direction, "reverse")
    
    def test_create_breakage(self):
        random.seed(1)
        self.EXL1.create_breakage()
        self.assertEqual(len(self.EXL1.exon_list), 2)
        self.assertEqual(self.EXL1.exon_list[1].get_start(), 800)
        self.assertEqual(self.EXL1.exon_list[1].get_end(), 1000)
        self.assertEqual(self.EXL1.junction, 999)
        self.setUp()
        random.seed(1)
        self.EXL1.create_breakage(1.0)
        self.assertEqual(len(self.EXL1.exon_list), 3)
        self.assertEqual(self.EXL1.exon_list[2].get_start(), 1500)
        self.assertEqual(self.EXL1.exon_list[2].get_end(), 1882)
        self.assertEqual(self.EXL1.junction, 1881)
        self.setUp()
        
    def test_get_possible_junction_exons(self):
        # if mid_exon_breaks allowed, all eons are possible, other wise the 
        # last or first exon are not allowed, depending on strand and whether
        # or not the strand is donor or acceptor
        self.assertEqual(self.EXL1.get_possible_junction_exons(False, 1, 1), 
                         ([0,1], False))
        self.assertEqual(self.EXL1.get_possible_junction_exons(True, 1, 1), 
                         ([0,1,2], True))
        self.assertEqual(self.EXL3.get_possible_junction_exons(False, 1, 1), 
                         ([1,2], False))
        self.assertEqual(self.EXL3.get_possible_junction_exons(True, 1, 1), 
                         ([0,1,2], True))
        # if no exons meet the brekable criteria, then no mid exon breakages 
        # allowed
        self.assertEqual(self.EXL1.get_possible_junction_exons(False, 1000, 1), 
                         ([0,1], False))
        self.assertEqual(self.EXL1.get_possible_junction_exons(True, 1000, 1), 
                         ([0,1], False))
        self.assertEqual(self.EXL3.get_possible_junction_exons(False, 1000, 1), 
                         ([1,2], False))
        self.assertEqual(self.EXL3.get_possible_junction_exons(True, 1000, 1), 
                         ([1,2], False))               
        
    
    def test_create_mid_exon_breakage(self):
        # exons are cut to the right of the junction exon, and the junction is
        # a random base in that exon   
        random.seed(1)
        self.EXL1.create_mid_exon_breakage(1, 1, 1)
        self.assertEqual(len(self.EXL1.exon_list), 2)
        self.assertEqual(self.EXL1.exon_list[1].get_start(), 800)
        self.assertEqual(self.EXL1.exon_list[1].get_end(), 827)
        self.assertEqual(self.EXL1.junction, 826)
        self.setUp()
        random.seed(1)
        self.EXL1.create_mid_exon_breakage(0, 1, 1)
        self.assertEqual(len(self.EXL1.exon_list), 1)
        self.assertEqual(self.EXL1.exon_list[0].get_start(), 100)
        self.assertEqual(self.EXL1.exon_list[0].get_end(), 154)
        self.assertEqual(self.EXL1.junction, 153)
        self.setUp()
        # exons are cut to the left of the junction exon, and the junction is
        # a random base in that exon
        random.seed(1)
        self.EXL3.create_mid_exon_breakage(1, 1, 1)
        self.assertEqual(len(self.EXL3.exon_list), 2)
        self.assertEqual(self.EXL3.exon_list[0].get_start(), 827)
        self.assertEqual(self.EXL3.exon_list[0].get_end(), 1000)
        self.assertEqual(self.EXL3.junction, 826)
        self.setUp()
        random.seed(1)
        self.EXL3.create_mid_exon_breakage(2, 1, 1)
        self.assertEqual(len(self.EXL3.exon_list), 1)
        self.assertEqual(self.EXL3.exon_list[0].get_start(), 1568)
        self.assertEqual(self.EXL3.exon_list[0].get_end(), 2000)
        self.assertEqual(self.EXL3.junction, 1567)
        self.setUp()
        
        
    def test_create_exon_breakage(self):
        # exons are cut to the right of the junction exon, and the junction is
        # the end of that exon
        self.EXL1.create_exon_breakage(1)
        self.assertEqual(len(self.EXL1.exon_list), 2)
        self.assertEqual(self.EXL1.exon_list[1].get_start(), 800)
        self.assertEqual(self.EXL1.junction, 999)
        self.setUp()
        self.EXL1.create_exon_breakage(0)
        self.assertEqual(len(self.EXL1.exon_list), 1)
        self.assertEqual(self.EXL1.exon_list[0].get_start(), 100)
        self.assertEqual(self.EXL1.junction, 499)
        self.setUp()
        # exons are cut to the left of the junction exon, and the junction is
        # the start of that exon
        self.EXL3.create_exon_breakage(1)
        self.assertEqual(len(self.EXL3.exon_list), 2)
        self.assertEqual(self.EXL3.exon_list[0].get_start(), 800)
        self.assertEqual(self.EXL3.junction, 799)
        self.setUp()
        self.EXL3.create_exon_breakage(2)
        self.assertEqual(len(self.EXL3.exon_list), 1)
        self.assertEqual(self.EXL3.exon_list[0].get_start(), 1500)
        self.assertEqual(self.EXL3.junction, 1499)
        self.setUp()
    
    def test_slice_exons_at_fusion_exon(self):
        # these exon lists remian unchanged
        self.EXL1.slice_exons_at_fusion_exon(2)
        self.assertEqual(len(self.EXL1.exon_list), len(self.exons1))
        self.EXL3.slice_exons_at_fusion_exon(0)
        self.assertEqual(len(self.EXL3.exon_list), len(self.exons1)) 
        # these are cut down to one exon
        self.EXL1.slice_exons_at_fusion_exon(0)
        self.assertEqual(len(self.EXL1.exon_list), 1)
        self.assertEqual(self.EXL1.exon_list[0].get_start(), 100)
        self.EXL3.slice_exons_at_fusion_exon(2)
        self.assertEqual(len(self.EXL3.exon_list), 1) 
        self.assertEqual(self.EXL3.exon_list[0].get_start(), 1500)
        self.setUp()
    
if __name__ == '__main__':
    unittest.main()