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
import requests

parser = argparse.ArgumentParser(description='Converts FUSIM output file to BEDPE format.')
parser.add_argument('fusim', help='Fusim txt output file.')
parser.add_argument('--bedpe', required=False, help='Output BEDPE', action='store_true')
parser.add_argument('--bed12', required=False, help='Output BED12', action='store_true')
args = parser.parse_args()


def getHGNC(ENSG):
	r = requests.get("http://grch37.rest.ensembl.org//xrefs/id/"+ENSG+"?external_db=HGNC", headers={ "Content-Type" : "application/json"})
	if not r.ok:
		r.raise_for_status()
 	decoded = r.json()[0]
 	return(decoded['display_id'])


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
				(nameA, nameB) = geneA[0].split('-')
				if geneB[4] == "+":
					geneBpos = min(geneB[8].split(','))
				else: 
					geneBpos = min(geneB[9].split(','))
				
				# subtract 1 from upstream position to be consistent with BED 0-based numbering	
				upstreamA = int(geneApos) - 1
				upstreamB = int(geneBpos) - 1
				bedpe = '\t'.join([geneA[3], str(upstreamA), geneApos, geneB[3], str(upstreamB), geneBpos, geneA[0], '0', geneA[4], geneB[4], getHGNC(nameA), getHGNC(nameB)])
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