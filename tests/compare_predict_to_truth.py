#! /usr/bin/env python

import sys

truthFusions = dict()
with open(sys.argv[1], 'r') as truth:
	for line in truth:
		if line.startswith('ENS'):
			vals = line.strip().split()
			truthFusions[str(vals[0])] = 0
			print vals[0]


print 'Number of fusions in truth: %s' % len(truthFusions)


with open(sys.argv[2], 'r') as predicted:
	for line in predicted:
		if not line.startswith('#'):
			vals = line.strip().split()
			gene1 = vals[4].split('^')[1]
			gene2 = vals[6].split('^')[1]
			fusedENSG = '-'.join([gene1, gene2])
			fusedENSGrev = '-'.join([gene2, gene1])
#			print '%s %s' % (fusedENSG, fusedENSGrev)
			if fusedENSG in truthFusions:
				truthFusions[fusedENSG] += 1
			elif fusedENSGrev in truthFusions:
				truthFusions[fusedENSGrev] += 1
				

nonZero = 0
for key,val in truthFusions.iteritems():
	if val > 0:
		print '%s: %s' % (key, val)
		nonZero += 1
		
print 'nonzero %s' % nonZero
