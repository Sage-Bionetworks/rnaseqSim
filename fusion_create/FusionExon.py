# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 14:03:41 2017

@author: aelamb
"""
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import random

class FusionExon(object):
    
    def __init__(self, seq_feature):
        self.seq_feature = seq_feature
    
    def get_seq_feature(self):
        return(self.seq_feature)
    
    def get_start(self):
        return(self.seq_feature.location.start)
    
    def get_end(self):
        return(self.seq_feature.location.end)
    
    def get_strand(self):
        return(self.seq_feature.location.strand)
    
    def set_start(self, loc):
        self.seq_feature.location = FeatureLocation(
            ExactPosition(loc), 
            self.get_end(),
            self.get_strand())
    
    def set_end(self, loc):
        self.seq_feature.location = FeatureLocation(
            self.get_start(),
            ExactPosition(loc), 
            self.get_strand())
    
    def __repr__(self):
        return(repr(self.seq_feature))
        
            
    def is_breakable(self, exon_min_size, exon_min_cleaved):
        exon_size = abs(self.get_start() - self.get_end())
        min_exon_left = exon_size - exon_min_cleaved
        return(min_exon_left >= exon_min_size)
    
    def break_exon_at_junction(self, exon_min_size, exon_min_cleaved, 
                               direction):
        if direction == "forward":
            first_possible_base = self.get_start() + exon_min_size
            last_possible_base = self.get_end() - exon_min_cleaved
        else:
            first_possible_base = self.get_start() + exon_min_cleaved
            last_possible_base = self.get_end() - exon_min_size
        junction_pos = random.randint(first_possible_base, last_possible_base)
        if direction == "forward":
            self.set_end(junction_pos)
        else:
            self.set_start(junction_pos)
        return(junction_pos)
    
        
