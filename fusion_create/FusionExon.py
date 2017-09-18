# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 14:03:41 2017

@author: aelamb
"""

from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import random

class FusionExon(SeqFeature):
    
    def __init__(self, seq_feature):
        self.location = seq_feature.location
        self.type = seq_feature.type
        self.id = seq_feature.id
        self.qualifiers = seq_feature.qualifiers
    
    def get_start(self):
        return(self.location.start)
    
    def get_end(self):
        return(self.location.end)
    
    def set_start(self, loc):
        self.location = FeatureLocation(
            ExactPosition(loc), 
            self.location.end,
            self.location.strand)
    
    def set_end(self, loc):
        self.location = FeatureLocation(
            self.location.start,
            ExactPosition(loc), 
            self.location.strand)
            
    def is_breakable(self, exon_min_size, exon_min_cleaved):
        exon_size = abs(self.location.start - self.location.end)
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
    
        
