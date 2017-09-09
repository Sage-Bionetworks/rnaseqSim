# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 07:34:53 2017

@author: aelamb
"""
#from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import FusionEventFunctions as FEV

class FusionEvent(dict):
    
    def add_mid_exon_fusion_breakpoints(
        self, 
        event_prob = 1.0, 
        two_break_prob = 0.5,
        left_break_prob = 0.5,
        min_bases_removed = 1,
        min_exon_size = 1):
            
        for exon_list in ['donorExons', 'acceptorExons']:
            self.shorten_exons(self[exon_list], event_prob, two_break_prob,
                left_break_prob, min_bases_removed, min_exon_size)
        
    @staticmethod
    def shorten_exons(exon_list, event_prob, two_break_prob, left_break_prob,
        min_bases_removed, min_exon_size):
                          
        for exon in exon_list:
            start, end = FEV.shorten_exon(exon, event_prob, two_break_prob, 
                left_break_prob, min_bases_removed, min_exon_size)
            exon.location = FeatureLocation(ExactPosition(start), 
                ExactPosition(end), strand = exon.location.strand)

                
    
    def add_sense_antisense_fusions(self):
        pass
    
    def add_exon_duplications_and_deletions(self):
        pass
    
    def add_fusion_events_in_UTR(self):
        pass
