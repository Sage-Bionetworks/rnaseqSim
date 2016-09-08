#! /usr/bin/env python
# KKD for Sage Bionetworks
# 4 Aug 2016

# sys.argv[1] = fusim TXT output
# sys.argv[2] = refflat to check against
# sys.argv[3] = 1-based position of transStart in refflat


import sys
import subprocess
import re

with open(sys.argv[1],'r') as fusions:
	for line in fusions:
		if not line.startswith('ENS'): continue
		vals = line.strip().split()
#		print '%s' % vals[2]
		refflatLine = subprocess.check_output(' '.join(['grep', vals[2], sys.argv[2]]), shell = True).strip()
		refflatSplit = refflatLine.split()
		transStart = refflatSplit[(int(sys.argv[3])-1)]
		transEnd = refflatSplit[int(sys.argv[3])]
		
		exonStarts = int(sys.argv[3])+5
		matched = 0
		for item in vals[exonStarts].split(','):
			if re.search(item,refflatLine) is not None:
				matched += 1
				break
		for item in vals[(exonStarts+1)].split(','):
			if re.search(item,refflatLine) is not None:
				matched += 1
				break
		if matched == 0:
			print 'no-match: %s' % vals[2]
