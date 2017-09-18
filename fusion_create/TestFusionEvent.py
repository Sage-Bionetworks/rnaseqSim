# -*- coding: utf-8 -*-

from FusionEvent import FusionEvent
from FusionExon import FusionExon
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import unittest
import random

class TestFusionEvent(unittest.TestCase):
    
    def setUp(self):
        random.seed(1)
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
        self.FE1 = FusionEvent(self.sfs1, self.sfs2, "+", "-")
        self.FE2 = FusionEvent(self.sfs1, self.sfs2, "-", "+")
    
    def test_getters(self):
        self.assertEqual(self.FE1.get_donor_exons().exon_list[0].get_start(), 
                         100)
        self.assertEqual(self.FE1.get_donor_exons().exon_list[2].get_end(), 
                         2000)
        self.assertEqual(self.FE1.get_acceptor_exons().exon_list[0].get_start(), 
                         10000)
        self.assertEqual(self.FE1.get_acceptor_exons().exon_list[2].get_end(), 
                         25000)
        self.assertEqual(self.FE1.get_donor_junction(), None) 
        self.assertEqual(self.FE1.get_acceptor_junction(), None) 
        
    def test_create_breakages1(self):
        # no mid-exon breakages
        # donor strand  looses its last exon since it is on the + strand
        self.FE1.create_breakages()
        self.assertEqual(len(self.FE1.get_donor_exons().exon_list), 2)
        self.assertEqual(self.FE1.get_donor_exons().exon_list[0].get_start(), 
                         100)
        self.assertEqual(self.FE1.get_donor_exons().exon_list[0].get_end(), 
                         500)
        self.assertEqual(self.FE1.get_donor_exons().exon_list[1].get_start(), 
                         800)
        self.assertEqual(self.FE1.get_donor_exons().exon_list[1].get_end(), 
                         1000)
        # acceptor strand  looses its first exon since it is on the - strand          
        self.assertEqual(len(self.FE1.get_acceptor_exons().exon_list), 2)
        self.assertEqual(self.FE1.get_acceptor_exons().exon_list[0].get_start(), 
                         17000)
        self.assertEqual(self.FE1.get_acceptor_exons().exon_list[0].get_end(), 
                         20000)
        self.assertEqual(self.FE1.get_acceptor_exons().exon_list[1].get_start(), 
                         22000)
        self.assertEqual(self.FE1.get_acceptor_exons().exon_list[1].get_end(), 
                         25000)
        
        # junctions are -1 from the end/start of the last/first exons
        self.assertEqual(self.FE1.get_donor_junction(), 999) 
        self.assertEqual(self.FE1.get_acceptor_junction(), 16999)
    
    def test_create_breakages2(self):
        # mid-exon breakages
        # donor strand  looses its first and last exon, and the start of its 
        # only exon breaks since it is on the - strand
        self.FE2.create_breakages(1.0)
        self.assertEqual(len(self.FE2.get_donor_exons().exon_list), 1)
        self.assertEqual(self.FE2.get_donor_exons().exon_list[0].get_start(), 
                         1882)
        self.assertEqual(self.FE2.get_donor_exons().exon_list[0].get_end(), 
                         2000)
        # acceptor strand  looses its first exon and break happens on  end of 
        # last exon since it is on the - strand
        self.assertEqual(len(self.FE2.get_acceptor_exons().exon_list), 2)
        self.assertEqual(self.FE2.get_acceptor_exons().exon_list[0].get_start(), 
                         10000)
        self.assertEqual(self.FE2.get_acceptor_exons().exon_list[0].get_end(), 
                         15000)
        self.assertEqual(self.FE2.get_acceptor_exons().exon_list[1].get_start(), 
                         17000)
        self.assertEqual(self.FE2.get_acceptor_exons().exon_list[1].get_end(), 
                         18349)
        
        # junctions are -1 from the end/start of the last/first exons
        self.assertEqual(self.FE2.get_donor_junction(), 1881) 
        self.assertEqual(self.FE2.get_acceptor_junction(), 18348) 
        self.setUp()
    

        
        
        
                          


    
#    def test_create_breakage(self):
#        # both strands are "-", donor gets cut off on left, acceptor on right
#        random.seed(1)
#        res1 = self.FE1.create_breakage("donor")
#        self.assertEqual(res1[0], 799)
#        self.assertEqual(res1[1], [self.exon2, self.exon3])
#        random.seed(1)
#        res2 = self.FE1.create_breakage("acceptor")
#        self.assertEqual(res2[0], 14999)
#        self.assertEqual(res2[1], [self.exon4])
#        # both strands are "+", donor gets cut off on right, acceptor on left
#        random.seed(1)
#        res3 = self.FE2.create_breakage("donor")
#        self.assertEqual(res3[0], 499)
#        self.assertEqual(res3[1], [self.exon1])
#        random.seed(4)
#        res4 = self.FE2.create_breakage("acceptor")
#        self.assertEqual(res4[0], 16999)
#        self.assertEqual(res4[1], [self.exon5, self.exon6])
#        
#        # mid exon breakages
#        random.seed(2)
#        res5 = self.FE1.create_breakage("donor", 
#                                        mid_exon_fusions = True, 
#                                        mid_exon_prob = 1.0)
#        self.assertEqual(res5[0], 1528)
#        self.assertEqual(res5[1][0].location.start, ExactPosition(1529))
#        self.assertEqual(res5[1][0].location.end, ExactPosition(2000))
#        random.seed(2)
#        res6 = self.FE1.create_breakage("acceptor", 
#                                        mid_exon_fusions = True,
#                                        mid_exon_prob = 1.0)
#        self.assertEqual(res6[0], 22169)
#        self.assertEqual(res6[1][0].location.start, ExactPosition(10000))
#        self.assertEqual(res6[1][0].location.end, ExactPosition(15000))
#        self.assertEqual(res6[1][1].location.start, ExactPosition(17000))
#        self.assertEqual(res6[1][1].location.end, ExactPosition(20000))                 
#        self.assertEqual(res6[1][2].location.start, ExactPosition(22000))
#        self.assertEqual(res6[1][2].location.end, ExactPosition(22170))                             

if __name__ == '__main__':
    unittest.main()