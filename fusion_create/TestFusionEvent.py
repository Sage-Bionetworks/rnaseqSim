# -*- coding: utf-8 -*-

from FusionEvent import FusionEvent
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import unittest

class TestFusionEvent(unittest.TestCase):
    
    def setUp(self):
        self.exon1 = SeqFeature(FeatureLocation(ExactPosition(95451328), 
                                          ExactPosition(95451419), 
                                          strand = -1),
                          type = 'exon', id = 'exon_770035')
        self.exon2 = SeqFeature(FeatureLocation(ExactPosition(131618), 
                                          ExactPosition(131645), 
                                          strand = -1), 
                          type = 'exon', id = 'exon_956247')
                          
        self.FE = FusionEvent({
              'dJunction': 95451328, 
              'aJunction': 169339, 
              'donorExons': [self.exon1], 
              'acceptorExons': [self.exon2]})
               
    def test_add_mid_exon_fusion_breakpoints1(self):
        # event_prob is 0, so nothing changes
        self.FE.add_mid_exon_fusion_breakpoints(event_prob = 0.0)
        self.assertEqual(self.FE['donorExons'][0].location.start, 95451328)
        self.assertEqual(self.FE['donorExons'][0].location.end, 95451419)
        self.assertEqual(self.FE['acceptorExons'][0].location.start, 131618)
        self.assertEqual(self.FE['acceptorExons'][0].location.end, 131645)
        self.setUp()
        # two_break_prob is 1.0, so all positions will change
        
    def test_add_mid_exon_fusion_breakpoints2(self):
        self.FE.add_mid_exon_fusion_breakpoints(event_prob = 1.0, 
                                                two_break_prob = 1.0)
        self.assertNotEqual(self.FE['donorExons'][0].location.start, 95451328)
        self.assertNotEqual(self.FE['donorExons'][0].location.end, 95451419)
        self.assertNotEqual(self.FE['acceptorExons'][0].location.start, 131618)
        self.assertNotEqual(self.FE['acceptorExons'][0].location.end, 131645)
        self.setUp()
    
    def test_add_mid_exon_fusion_breakpoints3(self):
        # two_break_prob is 0.0, left_break_prob is 1.0, so all starts will 
        # change
        self.FE.add_mid_exon_fusion_breakpoints(event_prob = 1.0, 
                                                two_break_prob = 0.0,
                                                left_break_prob = 1.0)
        self.assertNotEqual(self.FE['donorExons'][0].location.start, 95451328)
        self.assertEqual(self.FE['donorExons'][0].location.end, 95451419)
        self.assertNotEqual(self.FE['acceptorExons'][0].location.start, 131618)
        self.assertEqual(self.FE['acceptorExons'][0].location.end, 131645)
        self.setUp()
        
    def test_add_mid_exon_fusion_breakpoints4(self):
        # two_break_prob is 0.0, left_break_prob is 0.0, so all ends will 
        # change
        self.FE.add_mid_exon_fusion_breakpoints(event_prob = 1.0, 
                                                two_break_prob = 0.0,
                                                left_break_prob = 0.0)
        self.assertEqual(self.FE['donorExons'][0].location.start, 95451328)
        self.assertNotEqual(self.FE['donorExons'][0].location.end, 95451419)
        self.assertEqual(self.FE['acceptorExons'][0].location.start, 131618)
        self.assertNotEqual(self.FE['acceptorExons'][0].location.end, 131645)
        self.setUp()
        
  
 

if __name__ == '__main__':
    unittest.main()


