#! /usr/bin/env python
# KKD for Sage Bionetworks
# Mar 21, 2016

#fusionGene
#geneName
#name
#chrom
#strand
#exonCount
#exonBases
#exonIndexes
#exonStarts
#exonEnds
#fusionType
#fusionOptions

import argparse

parser = argparse.ArgumentParser(description='Converts FUSIM output file to BEDPE format.')
parser.add_argument('fusim', help='Fusim txt output file.')
parser.add_argument('--bedpe', required=False, help='Output BEDPE', action='store_true')
parser.add_argument('--bed12', required=False, help='Output BED12', action='store_true')
args = parser.parse_args()

with open(args.fusim, 'r') as fusim:
	geneA = None
	for line in fusim:

		if args.bedpe is True:
#			print line
			if not line.startswith('ENS'): continue
			if geneA is None:
				geneA = line.strip().split()
				if geneA[4] == "+":
					geneApos = max(geneA[9].split(','))
				else:
					geneApos = max(geneA[8].split(','))
#				print '%s' % geneApos
			else:
				geneB = line.strip().split()
				if geneB[4] == "+":
					geneBpos = min(geneB[8].split(','))
				else: 
					geneBpos = min(geneB[9].split(','))
				bedpe = '\t'.join([geneA[3], geneApos, geneApos, geneB[3], geneBpos, geneBpos, geneA[0], '0', geneA[4], geneB[4]])
				print '%s' % bedpe
				geneA = None
				geneB = None
			
			
	#	if args.bed12 is True:		
	# 		if not line.startswith('ENSG'): continue
	# 		vals = line.strip().split
	# 		start = min(vals[8].split(','))
	# 		end = max(vals[9].split(','))
	# 		bed12 = '\t'.join([vals[3], start, end, vals[0], 0, valsA[4], ])
	# 		print '%s' % bedpe
	# 		valsA = ''
	# 		valsB = ''