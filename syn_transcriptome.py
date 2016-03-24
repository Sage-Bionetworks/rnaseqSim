#! /usr/bin/env python
# Original version written by Yin Hu in 2015 for Sage Bionetworks

# This script synthesizes the transcriptome
# - sub-samples from the complete set of annotated transcripts
# - calculates the transcript lengths

import csv, random

expressprob = 0.3

with open('Fusim//fusions.fasta', 'r') as fafile, open('translen_fusion.txt', 'w') as lenfile:
	reader = csv.reader(fafile, delimiter="\t")
	for l1 in reader:
		l2 = reader.next()
		id = str(l1[0]).split()[0]
		#id = id[id.rfind('|')+1:]
		lenfile.write(id + "\t" + str(len(str(l2[0]))) + "\n")

with open('Ref/gencode.v19.pc_transcripts.fa', 'r') as fafile, open('ref.fa', 'w') as samplefile, open('translen_ref.txt', 'w') as lenfile:
	reader = csv.reader(fafile, delimiter="\t")
	for l1 in reader:
		l2 = reader.next()
		
		if random.random() > expressprob:
			continue
			
		id = str(l1[0]).split()[0]
		samplefile.write(id + "\n")
		samplefile.write(str(l2[0]) + "\n")
		lenfile.write(id + "\t" + str(len(str(l2[0]))) + "\n")



