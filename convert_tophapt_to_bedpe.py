#! /usr/bin/env python
# convert tophat fusions.out to BEDPE


import sys

with open(sys.argv[1], 'r') as tophat:
	for line in tophat:
		vals = line.strip().split()
		chroms = vals[0].split('-')
		name = '-'.join([vals[0], vals[1], vals[2]])
		start1 = int(vals[1]) - 1
		start2 = int(vals[2]) -1
		if vals[3] == 'ff': 
			strand1 = '+'
			strand2 = '+'
		elif vals[3] == 'fr': 
			strand1 = '+'
			strand2 = '-'
		elif vals[3] == 'rf': 
			strand1 = '-'
			strand2 = '+'
		elif vals[3] == 'rr': 
			strand1 = '-'
			strand2 = '-'
		else:
			print 'strange strandedness'
		bedpe = '\t'.join([chroms[0], str(start1), vals[1], chroms[1], str(start2), vals[2], name, '0', strand1, strand2 ])
		print '%s' % bedpe