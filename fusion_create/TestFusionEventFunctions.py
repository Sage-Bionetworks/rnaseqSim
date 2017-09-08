# -*- coding: utf-8 -*-
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import FusionEventFunctions as FEV
import unittest

class TestFusionEventFunctions(unittest.TestCase):
    
    def setUp(self):
        self.exon1 = SeqFeature(FeatureLocation(ExactPosition(95451328), 
                                                ExactPosition(95451419), 
                                                strand = -1),
                          type = 'exon', id = 'exon_770035')
                          
    def test_shorten_exon(self):
        # positions not changed due to probs
        self.assertEqual(FEV.shorten_exon(
            self.exon1, 0.0, 0.0, 0.0, 1, 80)[0], 95451328)
        self.assertEqual(FEV.shorten_exon(
            self.exon1, 0.0, 0.0, 0.0, 1, 80)[1], 95451419)
        # positions not changed due to min exon size
        self.assertEqual(FEV.shorten_exon(
            self.exon1, 1.0, 0.0, 0.0, 1, 91)[0], 95451328)
        self.assertEqual(FEV.shorten_exon(
            self.exon1, 1.0, 0.0, 0.0, 1, 91)[1], 95451419)
        # left break occurs
        self.assertNotEqual(FEV.shorten_exon(
            self.exon1, 1.0, 0.0, 1.0, 1, 90)[0], 95451328)
        self.assertEqual(FEV.shorten_exon(
            self.exon1, 1.0, 0.0, 1.0, 1, 90)[1], 95451419)
        # right break occurs
        self.assertEqual(FEV.shorten_exon(
            self.exon1, 1.0, 0.0, 0.0, 1, 90)[0], 95451328)
        self.assertNotEqual(FEV.shorten_exon(
            self.exon1, 1.0, 0.0, 0.0, 1, 90)[1], 95451419)
        # right break occurs beacuse min exon size two large for two breaks
        self.assertEqual(FEV.shorten_exon(
            self.exon1, 1.0, 1.0, 0.0, 1, 90)[0], 95451328)
        self.assertNotEqual(FEV.shorten_exon(
            self.exon1, 1.0, 1.0, 0.0, 1, 90)[1], 95451419)
        # two breaks occur
        self.assertNotEqual(FEV.shorten_exon(
            self.exon1, 1.0, 1.0, 0.0, 1, 89)[0], 95451328)
        self.assertNotEqual(FEV.shorten_exon(
            self.exon1, 1.0, 1.0, 0.0, 1, 89)[1], 95451419)

        
                          
    def test_determine_exon_event(self):
        # max exon size is too large, so probs don't matter
        self.assertEqual(FEV.determine_exon_event(
            self.exon1, 0.0, 0.0, 0.0, 1, 91), 'none')
        self.assertEqual(FEV.determine_exon_event(
            self.exon1, 1.0, 1.0, 1.0, 1, 91), 'none')
        # max exon size is too large for two break events
        self.assertEqual(FEV.determine_exon_event(
            self.exon1, 0.0, 0.0, 0.0, 1, 90), 'none')
        self.assertEqual(FEV.determine_exon_event(
            self.exon1, 1.0, 0.0, 0.0, 1, 90), 'right')
        self.assertEqual(FEV.determine_exon_event(
            self.exon1, 1.0, 0.0, 1.0, 1, 90), 'left')
        self.assertEqual(FEV.determine_exon_event(
            self.exon1, 1.0, 1.0, 1.0, 1, 90), 'left')
        # one and two break events can happen
        self.assertEqual(FEV.determine_exon_event(
            self.exon1, 0.0, 0.0, 0.0, 1, 89), 'none')
        self.assertEqual(FEV.determine_exon_event(
            self.exon1, 1.0, 0.0, 0.0, 1, 89), 'right')
        self.assertEqual(FEV.determine_exon_event(
            self.exon1, 1.0, 0.0, 1.0, 1, 89), 'left')
        self.assertEqual(FEV.determine_exon_event(
            self.exon1, 1.0, 1.0, 1.0, 1, 89), 'both')
            
        
    def test_determine_exon_max_breaks(self):
        self.assertEqual(FEV.determine_exon_max_breaks(
            self.exon1, 1, 91), 0)
        self.assertEqual(FEV.determine_exon_max_breaks(
            self.exon1, 1, 90), 1)
        self.assertEqual(FEV.determine_exon_max_breaks(
            self.exon1, 1, 89), 2)
        self.assertEqual(FEV.determine_exon_max_breaks(
            self.exon1, 21, 71), 0)
        self.assertEqual(FEV.determine_exon_max_breaks(
            self.exon1, 20, 71), 1)
        self.assertEqual(FEV.determine_exon_max_breaks(
            self.exon1, 10, 71), 2)
        
    def test_shorten_exon_on_start(self):
        self.assertNotEqual(FEV.shorten_exon_on_start(
            self.exon1, 1, 1)[0], 95451328)
        self.assertEqual(FEV.shorten_exon_on_start(
            self.exon1, 1, 1)[1], 95451419)
    
    def test_shorten_exon_on_end(self):
        self.assertEqual(FEV.shorten_exon_on_end(
            self.exon1, 1, 1)[0], 95451328)
        self.assertNotEqual(FEV.shorten_exon_on_end(
            self.exon1, 1, 1)[1], 95451419)
    
    def test_shorten_exon_on_both_sides(self):
        self.assertNotEqual(FEV.shorten_exon_on_both_sides(
            self.exon1, 1, 1)[0], 95451328)
        self.assertNotEqual(FEV.shorten_exon_on_both_sides(
            self.exon1, 1, 1)[1], 95451419)
        self.assertEqual(FEV.shorten_exon_on_both_sides(
            self.exon1, 1, 89)[0], 95451329)
        self.assertEqual(FEV.shorten_exon_on_both_sides(
            self.exon1, 1, 89)[1], 95451418)
        

if __name__ == '__main__':
    unittest.main()