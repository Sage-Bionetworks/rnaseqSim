#! /usr/bin/env python3
# KKD for Sage Bionetworks
# 25 Aug. 2016


import sys

firstStart = 0
otherStart = 0
geneA = None
geneB = ''
lines = 0
with open(sys.argv[1], 'r') as fusim:
	for line in fusim:
		if not line.startswith('ENSG'): continue
		if geneA is None:
			geneA = line.strip().split()
			lines += 1
		else:
			geneB = line.strip().split()
			assert geneA[0] == geneB[0]
			if geneB[4] == "+":
				if geneB[7].startswith('0'):
					firstStart += 1
				else:
					otherStart += 1
			geneA = None
			lines += 1
print('firstStart: '+str(firstStart)+' otherstart: '+str(otherStart))
print('lines: '+str(lines))