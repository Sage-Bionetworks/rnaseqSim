#! /usr/bin/env python
# KKD for Sage Bionetworks
# Mar 11, 2016

import argparse

parser = argparse.ArgumentParser(description='Filters output of Fusim by transcript length')
parser.add_argument('fusim', help='Fusim FASTA output file.')
parser.add_argument('--length', required=False, help='Minimum length for fused transcripts.', default = 400, type=int)
parser.add_argument('--number', required=False, help='Number of transcripts needed that meet the length threshold.', default = 10, type=int)
args = parser.parse_args()



numberMeetingCriteria = 0
with open(args.fusim, 'r') as fusim:
	outfile = ''.join([args.fusim.rstrip('.fasta'), '_filtered', '.gtf'])
	with open(outfile, 'w') as gtf:
#		gtf.write('# GTF-formatted fusim output\n')
		for line in fusim:
			if numberMeetingCriteria >= args.number: break
			elif line.startswith('>'):
				fullName = line
			else:
				length = len(line.strip())
				if length > args.length:
					print fullName,line
					numberMeetingCriteria += 1
					
					seqname = fullName.strip().split()[0].lstrip('>ref|')
					field9 = ' '.join(['gene_id "'+seqname+'";', 'gene_biotype "fusion";'])
					geneGtfLine = '\t'.join([seqname, 'fusim', 'gene', '1', str(length), '.', '+', '.', field9])
					gtf.write(geneGtfLine+'\n')
					field9 = ' '.join(['gene_id "'+seqname+'";', 'transcript_id "'+seqname+'";', 'gene_biotype "fusion";'])
					transcriptGtfLine = '\t'.join([seqname, 'fusim', 'transcript', '1', str(length), '.', '+', '.', field9])
					gtf.write(transcriptGtfLine+'\n')
					field9 = ' '.join(['gene_id "'+seqname+'";', 'transcript_id "'+seqname+'";', 'exon_number "1";', 'gene_biotype "fusion";'])
					gtfLine = '\t'.join([seqname, 'fusim', 'exon', '1', str(length), '.', '+', '.', field9])
					gtf.write(gtfLine+'\n')
