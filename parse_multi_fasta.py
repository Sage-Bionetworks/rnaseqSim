#! /usr/bin/env python


import os, argparse

parser = argparse.ArgumentParser(description='Split a multi-fasta into individual fasta files.')
parser.add_argument('fasta', help='Input multi-fasta.')
args = parser.parse_args()


with open(args.fasta, 'r') as input:
	for line in input:
		if line.startswith('>'):
			vals = line.strip().split()
			if output:
				output.close()
			output = open(vals[0].lstrip('>'),'w')
		print output.write(line)
