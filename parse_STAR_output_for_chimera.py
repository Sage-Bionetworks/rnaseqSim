#! /usr/bin/env python
# KKD for Sage Bionetworks
# Dec. 1, 2015
# Parses star chimeric.out.junction output files for input format needed by 'chimera' r package


import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Submits a workflow job to a grid queue manager such as SGE.')
parser.add_argument('starFile', help='STAR output file \'*Chimeric.out.junction.\'')
parser.add_argument('--excludeMT', required=False, help='Exclude MT chimera from output?', action='store_true')
args = parser.parse_args()


breakpoints = dict()

with open(args.starFile) as starInput:
	for line in starInput:
		vals = line.strip().split()
		chimType = vals[6]
		if chimType == "-1": continue
		
		gene1Chr = vals[0]
		gene1Pos = vals[1]
		gene1Str = vals[2]
		gene2Chr = vals[3]
		gene2Pos = vals[4]
		gene2Str = vals[5]
		junction = "_".join([gene1Chr,gene1Pos,gene1Str,gene2Chr,gene2Pos,gene2Str])
		if args.excludeMT is True and (gene1Chr == "MT" or gene1Chr == "M" or gene2Chr == "MT" or gene2Chr == "M"): continue

		if junction not in breakpoints:
			breakpoints[junction] = [1,0]
		else:
			breakpoints[junction][0] += 1
			
starInput.close()			


sortedBreakpoints = dict()
for junction in breakpoints:
	vals = junction.split("_")
	if vals[0] not in sortedBreakpoints:
	# 0= gene1Pos, 1=gene1Str, 2=gene2Chr, 3=gene2Pos, 4=gene2Str
		sortedBreakpoints[vals[0]] = [(vals[1],vals[2],vals[3],vals[4],vals[5])]
	else: 
		sortedBreakpoints[vals[0]].append((vals[1],vals[2],vals[3],vals[4],vals[5]))
#print 'Number of chromsomes in sorted dictionary: %d' % len(sortedBreakpoints)
	
for chr in sortedBreakpoints:
	sortedBreakpoints[chr].sort()
#	print 'Size of chr %s list: %d' % (chr, len(sortedBreakpoints[chr]))


obsCount = 0
with open(args.starFile) as starInput:
	for line in starInput:
		vals = line.strip().split()
		chimType = vals[6]
		if chimType != "-1": continue
		
		gene1Chr = vals[0]
		gene1Pos = vals[1]
		gene2Chr = vals[3]
		gene2Pos = vals[4]
		
		if gene1Chr not in sortedBreakpoints:
			obsCount += 1
			continue
		for junction in sortedBreakpoints[gene1Chr]:
			if gene1Pos > junction[0]: break
			# currently counting encompassing reads regardless of strand matching
			elif gene1Pos < junction[0] and gene2Chr == junction[2] and gene2Pos > junction[3]:
				breakpoints["_".join([gene1Chr,junction[0],junction[1],gene2Chr,junction[3],junction[4]])][1] += 1
starInput.close()			

print 'Number of chimeric encompassing reads that are only observed on scaffolds and only as encompassing: %d' % obsCount


for junction in breakpoints:
	vals = junction.split("_")
#	print junction, breakpoints[junction][0], breakpoints[junction][1]	
	print '%s\t%d\t%d' % ('\t'.join(map(str,vals)), breakpoints[junction][0], breakpoints[junction][1])