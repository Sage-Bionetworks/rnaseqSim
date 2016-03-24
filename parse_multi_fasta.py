#! /usr/bin/env python


import argparse

parser = argparse.ArgumentParser(description='Split a multi-fasta into individual fasta files.')
parser.add_argument('fasta', help='Input multi-fasta.')
args = parser.parse_args()


with open(args.fasta, 'r') as input:
	for line in input:
		if line.startswith('>'):
			vals = line.strip().split()
			try:
				output.close()
			except NameError:
				pass
			output = open('.'.join([vals[0].lstrip('>ref|'), 'fasta']),'w')
			output.write('>'+vals[0].lstrip('>ref|')+'\n')
		else:
			output.write(line)
