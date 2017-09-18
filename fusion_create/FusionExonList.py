# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:09:42 2017

@author: aelamb
"""

from FusionExon import FusionExon
import random

class FusionExonList(object):
    
    def __init__(self, lst, strand, category):
        
        self.exon_list = [FusionExon(exon) for exon in lst]
        self.set_direction(strand, category)
        self.junction = None
        
    def get_exon_list(self):
        return(self.exon_list)
    
    def get_junction(self):
        return(self.junction)
        
    def __len__(self):
        return(len(self.exon_list))
    
    def __getitem__(self, position):
        return(self.exon_list[position])
    
    def __repr__(self):
        return("FusionExonList(%r, %r, %r)" % 
            (self.exon_list, self.direction, self.junction))
        
    def set_direction(self, strand, category):
        if ((category == "donor" and strand == "+") or 
            (category == "acceptor" and strand == "-")):
            self.direction = "forward"
        else:
            self.direction = "reverse"
    
    def create_breakage(self,
                        mid_exon_prob = 0.0,
                        exon_min_size = 1, 
                        exon_min_cleaved = 1):
                            
        allow_mid_exon = random.random() <= mid_exon_prob
        junction_exons, do_mid_exon = self.get_possible_junction_exons(
            allow_mid_exon, exon_min_size, exon_min_cleaved)
        junction_exon_n = random.choice(junction_exons)
        if do_mid_exon:
            self.create_mid_exon_breakage(
                junction_exon_n, exon_min_size, exon_min_cleaved)
        else:
            self.create_exon_breakage(junction_exon_n)

    def get_possible_junction_exons(self, allow_mid_exon, exon_min_size, 
                                    exon_min_cleaved):
        idx_range = range(0, len(self))
        if allow_mid_exon:
            idx_lst = [idx for idx in idx_range if 
                self[idx].is_breakable(exon_min_size, exon_min_cleaved)]
            if len(idx_lst) == 0:
                allow_mid_exon = False
        if not allow_mid_exon:
            if self.direction == "forward":
                idx_lst = idx_range[:-1]
            else:
                idx_lst = idx_range[1:]
        return(idx_lst, allow_mid_exon)
    
    def create_mid_exon_breakage(self, junction_exon_n, exon_min_size, 
                                 exon_min_cleaved):
        junc_exon = self[junction_exon_n]
        self.slice_exons_at_fusion_exon(junction_exon_n)
        self.junction = junc_exon.break_exon_at_junction(
            exon_min_size, exon_min_cleaved, self.direction) - 1
    
    def create_exon_breakage(self, junction_exon_n):
        junc_exon = self[junction_exon_n]
        self.slice_exons_at_fusion_exon(junction_exon_n)
        if self.direction == "forward":
            junction = junc_exon.get_end()
        else:
            junction = junc_exon.get_start()
        self.junction = junction - 1 
        
    def slice_exons_at_fusion_exon(self, junction_exon_n):
        if self.direction == "forward":
            self.exon_list = self[:junction_exon_n + 1]
        else:
            self.exon_list = self[junction_exon_n:]         