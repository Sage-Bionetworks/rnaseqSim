#! /usr/bin/env python
# KKD for Sage Bionetworks
# June 13, 2016
# Convert star-fusion output to BEDPE

import sys

with open(sys.argv[1], 'r') as starOutput:
	for line in starOutput:
		if not line.startswith('#'):
			vals = line.strip().split()
			valsL = vals[5].split(':')
			valsR = vals[7].split(':')
			chrL = valsL[0]
			posL = valsL[1]
			geneL = vals[4].split('^')[1]
			strandL = valsL[2]
			chrR = valsR[0]
			posR = valsR[1]
			geneR = vals[6].split('^')[1]
			strandR = valsR[2]
			bedpe = '\t'.join([chrL, str(int(posL)-1), posL, chrR, str(int(posR)-1), posR,'-'.join([geneL, geneR]), '0', strandL, strandR])
			print bedpe