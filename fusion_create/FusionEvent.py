# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 07:34:53 2017

@author: aelamb
"""
from FusionExonList import FusionExonList

class FusionEvent(object):
    
    def __init__(self, donor_features, acceptor_features, donor_strand, 
                 acceptor_strand):
        self.donor_exons = FusionExonList(donor_features, donor_strand, 
                                          "donor")
        self.acceptor_exons = FusionExonList(acceptor_features, donor_strand, 
                                             "acceptor")
    
    def get_donor_exons(self):
        return(self.donor_exons)
    
    def get_acceptor_exons(self):
        return(self.acceptor_exons)
    
    def get_donor_junction(self):
        return(self.donor_exons.get_junction())
    
    def get_acceptor_junction(self):
        return(self.acceptor_exons.get_junction())
    
    def create_breakages(self, 
                         mid_exon_prob = 0.0,
                         min_exon_min_size = 1, 
                         min_exon_min_cleaved = 1):
        self.donor_exons.create_breakage(
            mid_exon_prob, min_exon_min_size, min_exon_min_cleaved)
        self.acceptor_exons.create_breakage(
            mid_exon_prob, min_exon_min_size, min_exon_min_cleaved)