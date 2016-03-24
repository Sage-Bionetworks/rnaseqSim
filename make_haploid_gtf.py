#! /usr/bin/env python
# Mar 7, 2016
# KKD for Sage Bionetworks


import sys
import argparse

parser = argparse.ArgumentParser(description='Adds suffix to chrom ids, gene ids, and transcript ids in input GTF to facilitate creation of diploid genome.')
parser.add_argument('GTF', help='GTF input file.')
parser.add_argument('suffix', help='Suffix label to add, e.g. "hap1"')
args = parser.parse_args()


with open(args.GTF, 'r') as gtf:
	for line in gtf:
		if not line.startswith('#'):
			main_vals = line.strip().split('\t')
			main_vals[0] = '-'.join([main_vals[0], args.suffix])
			if not len(main_vals) == 9:
<<<<<<< HEAD
				print '%s\n' % len(main_vals)
			else:
				attributes = main_vals.pop().strip().split(';')
#				attributes = main_vals[8].strip().split(';')
				for i in range(len(attributes)):
					if attributes[i].startswith('gene_id'):
						gene = attributes[i].split()
						gene[1] = '-'.join([gene[1].rstrip('"'), args.suffix+'"'])
						attributes[i] = ' '.join(gene)
					elif attributes[i].startswith(' transcript_id'):
						transcript = attributes[i].split()
						transcript[1] = '-'.join([transcript[1].rstrip('"'), args.suffix+'"'])
						attributes[i] = ' '.join(transcript)
			first_part = '\t'.join(main_vals)
			second_part = '; '.join(attributes)
			print '%s\t%s' % (first_part, second_part)
		else:
			print line.strip()
=======
				print line
			else:
				attributes = main_vals[8].strip().split(';')
				for i in range(len(attribues)):
					if attributes[i].startswith('gene_id'):
						gene = attributes[i].split('\s')
						gene[1] = '-'.join([gene[1].rstrip('"'), args.suffix+'"'])
						attributes[i] = ' '.join(gene)
					elif item.startswith('transcript_id'):
						transcript = item.split('\s')
						transcript[1] = '-'.join([transcript[1].rstrip('"'), args.suffix+'"'])
						attributes[i] = ' '.join(transcript)
			first_part = '\t'.join(main_vals)
			second_part = '; '.attributes.join(attributes)
			print '%s\t%s' % (first_part, second_part)
>>>>>>> 688fc3583f054923a9a56ffe6c009fdb0ffcbee3
