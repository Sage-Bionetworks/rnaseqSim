#! /usr/bin/env python


<<<<<<< HEAD
import argparse
=======
import os, argparse
>>>>>>> 688fc3583f054923a9a56ffe6c009fdb0ffcbee3

parser = argparse.ArgumentParser(description='Split a multi-fasta into individual fasta files.')
parser.add_argument('fasta', help='Input multi-fasta.')
args = parser.parse_args()


with open(args.fasta, 'r') as input:
	for line in input:
		if line.startswith('>'):
			vals = line.strip().split()
<<<<<<< HEAD
			try:
				output.close()
			except NameError:
				pass
			output = open('.'.join([vals[0].lstrip('>ref|'), 'fasta']),'w')
			output.write('>'+vals[0].lstrip('>ref|')+'\n')
		else:
			output.write(line)
=======
			if output:
				output.close()
			output = open(vals[0].lstrip('>'),'w')
		print output.write(line)
>>>>>>> 688fc3583f054923a9a56ffe6c009fdb0ffcbee3
