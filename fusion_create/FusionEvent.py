# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 07:34:53 2017

@author: aelamb
"""
import FusionEventFunctions as FEV
import random

class FusionEvent(object):
    
    def __init__(self, donor_exons, acceptor_exons, donor_strand, 
                 acceptor_strand):
        self.donor_exons = donor_exons
        self.acceptor_exons = acceptor_exons
        self.donor_strand = donor_strand
        self.acceptor_strand = acceptor_strand
        self.donor_junction = None
        self.acceptor_junction = None
    
    def get_donor_exons(self):
        return(self.donor_exons)
    
    def get_acceptor_exons(self):
        return(self.acceptor_exons)

    def get_donor_strand(self):
        return(self.donor_strand)
    
    def get_acceptor_strand(self):
        return(self.acceptor_strand)
    
    def get_donor_junction(self):
        return(self.donor_junction)
    
    def get_acceptor_junction(self):
        return(self.acceptor_junction)
    
    def create_breakages(self, mid_exon_fusions = False):
        self.donor_junction, self.donor_exons = self.create_breakage(
            "donor", mid_exon_fusions)
        self.acceptor_junction, self.acceptor_exons = self.create_breakage(
            "acceptor", mid_exon_fusions)
        
    def create_breakage(self, strand, mid_exon_fusions = False):
        if strand == "donor":
            exons = self.donor_exons
            strand_dir = self.donor_strand
        else:
            exons = self.acceptor_exons
            strand_dir = self.acceptor_strand
        direction = FEV.get_direction(strand, strand_dir)
        junction_range = FEV.get_junction_range(
            exons, direction, mid_exon_fusions)
        junction_exon_n = random.randint(*junction_range)
        if mid_exon_fusions:
            junction, exons = FEV.create_mid_exon_breakage(
                exons, junction_exon_n, direction)
        else:
            junction, exons = FEV.create_exon_breakage(
                exons, junction_exon_n, direction)
        junction -= 1
        return(junction, exons)

#    def __repr__(self):
#        pass
#        return(
#            self.get_donor_exons() +
#            self.get_acceptor_exons() +
#            self.get_donor_junction() +
#            self.get_acceptor_junction())
    
