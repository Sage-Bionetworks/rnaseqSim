# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 15:13:52 2017

@author: aelamb
"""
import random

def shorten_exon(exon, 
                 event_prob, 
                 two_break_prob, 
                 left_break_prob,
                 min_bases_removed,
                 min_exon_size):
                          
    event = determine_exon_event(exon, 
                                 event_prob, 
                                 two_break_prob, 
                                 left_break_prob,
                                 min_bases_removed,
                                 min_exon_size)     
    if event == "none":
        # exon will not be shortened
        start = exon.location.start
        end =  exon.location.end
    elif event == "both":
        # exon will be shortened on both sides
        start, end = shorten_exon_on_both_sides(
                exon, min_bases_removed, min_exon_size)
    elif event == "left":
        # exon will be shortened on left side
        start, end = shorten_exon_on_start(
                exon, min_bases_removed, min_exon_size)
    else:
        # exon will be shortened on right side
        start, end = shorten_exon_on_end(
                exon, min_bases_removed, min_exon_size)   
    return(start, end)

def determine_exon_event(exon, 
                         event_prob, 
                         two_break_prob, 
                         left_break_prob,
                         min_bases_removed,
                         min_exon_size):
    max_breaks = determine_exon_max_breaks(
        exon, min_bases_removed, min_exon_size)
    if random.random() > event_prob or max_breaks == 0:
        return("none")
    elif random.random() <= two_break_prob and max_breaks == 2:
        return("both")
    elif random.random() <= left_break_prob:
        return("left")
    else:
        return("right")

def determine_exon_max_breaks(exon, min_bases_removed, min_exon_size):
    exon_length = abs(exon.location.end - exon.location.start)
    if exon_length - (2 * min_bases_removed + min_exon_size) >= 0:
        return(2)
    elif exon_length - (min_bases_removed + min_exon_size) >= 0:
        return(1)
    else:
        return(0)

                
def shorten_exon_on_start(exon, min_bases_removed, min_exon_size):
    exon_break = random.randint(exon.location.start + min_bases_removed, 
                                exon.location.end - min_exon_size)
    return(exon_break, exon.location.end)

def shorten_exon_on_end(exon, min_bases_removed, min_exon_size):
    exon_break = random.randint(exon.location.start + min_exon_size, 
                                exon.location.end - min_bases_removed)
    return(exon.location.start, exon_break)
        
def shorten_exon_on_both_sides(exon, min_bases_removed, min_exon_size):
    break1_min = exon.location.start + min_bases_removed
    break1_max = exon.location.end - (min_exon_size + min_bases_removed)
    exon_break1 = random.randint(break1_min, break1_max)
    break2_min = exon_break1 + min_exon_size
    break2_max = exon.location.end  - min_bases_removed
    exon_break2 = random.randint(break2_min, break2_max)
    exon_breaks = [exon_break1, exon_break2]                       
    return(exon_breaks)