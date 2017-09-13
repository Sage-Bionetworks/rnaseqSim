# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 15:13:52 2017

@author: aelamb
"""
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
import random


def get_direction(strand, strand_dir):
    if ((strand == "donor" and strand_dir == "+") or 
        (strand == "acceptor" and strand_dir == "-")):
        direction = "forward"
    else:
        direction = "reverse"
    return(direction)

def get_junction_range(exons, direction, mid_exon_fusions):
     first_possible_exon_junction = 0
     last_possible_exon_junction = len(exons) - 1
     if not mid_exon_fusions:
         if direction == "forward":
             last_possible_exon_junction -= 1 
         else:
             first_possible_exon_junction += 1
     return(first_possible_exon_junction, last_possible_exon_junction)

def create_mid_exon_breakage(exons, junction_exon_n, direction):
    junc_exon = exons[junction_exon_n]
    first_possible_base = junc_exon.location.start + 1
    last_possible_base = junc_exon.location.end - 1
    junction = ExactPosition(random.randint(
        first_possible_base, last_possible_base))
    if direction == "forward":
        exons[junction_exon_n] = create_new_exon(junc_exon, "end", junction)
    else:
        exons[junction_exon_n] = create_new_exon(junc_exon, "start", junction)
    exons = slice_exons_at_fusion_exon(exons, junction_exon_n, direction)
    return(junction, exons)

def create_exon_breakage(exons, junction_exon_n, direction):
    if direction == "forward":
        junction = exons[junction_exon_n].location.end
    else:
        junction = exons[junction_exon_n].location.start
    exons = slice_exons_at_fusion_exon(exons, junction_exon_n, direction)
    return(junction, exons)
    
def slice_exons_at_fusion_exon(exons, junction_exon_n, direction):
    if direction == "forward":
        exons = exons[0:junction_exon_n + 1]
    else:
        exons = exons[junction_exon_n: len(exons)] 
    return(exons)

def create_new_exon(exon, loc_name, new_position):
    if loc_name == 'start':
        start = new_position
        end = exon.location.end
    else:
        start = exon.location.start
        end = new_position
    new_loc = FeatureLocation(ExactPosition(start), 
                              ExactPosition(end), 
                              exon.location.strand)                       
    return(SeqFeature(
        new_loc, type = 'exon', id = exon.id, qualifiers = exon.qualifiers))