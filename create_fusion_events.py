#! /usr/bin/env python

#from toil.job import Job
import subprocess
import argparse
import requests
import os
import re
from Bio import SeqIO
from Bio.Alphabet import IUPAC

    

referenceGenome = '/external-data/Genome/genomes/Hsapiens_Ensembl_GRCh37/primary_nomask/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
referenceGeneModels = '/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/Hsapiens_Ensembl_v75.genePred_refflat3'

########################
## Workflow functions
########################

def runFusim(numEvents, simName, memory="2G", cores=1, disk="1M"):
	'''Runs Fusim to generate fusion events.'''
	
	fusimJarPath = '/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Scripts/fusim-0.2.2/fusim.jar'
	
	numEventsToSimulate = int(numEvents)*5
	cmd = ''.join(['java -jar ', fusimJarPath, ' --gene-model=', referenceGeneModels, ' --fusions=', str(numEventsToSimulate), ' --reference=', referenceGenome, ' --fasta-output=', simName+'.fasta', ' --text-output=', simName+'fusions.txt', ' --auto-correct-orientation --cds-only --keep-exon-boundary'])
	print(cmd)
	subprocess.call(cmd, shell = True)


def filterFusionEvents(numEvents, minLength, simName, memory="100M", cores=1, disk="1M"):
	'''Reads in Fusim output files and writes out a GTF, FASTAs, and BEDPE of events meeting min length criteria.'''
	
	filteredFusimEvents = list()
	allFastaFilenames = '' # Needed for RSEM reference creation later

	# Get exon start/end sites from TXT output file
	moreAnnotations = parseFusimTxtOutput(simName+'fusions.txt')

	inputFusim = simName+'.fasta'				
	for record in SeqIO.parse(inputFusim, "fasta") :
		if len(filteredFusimEvents) == numEvents: break
		temp = record.id.split("|")[1] # Specific to FUSIM inputs
		record.id = temp
		record.annotations.update(parseFusimFastaDescription(record.description))
		record.annotations.update(moreAnnotations[record.id])

		# Filter out acceptor genes that start with exon 0 for + strand acceptors,
		# or donor genes that start (end) with 0 for - strand.
		if record.annotations['exonIndex2'].startswith('0') and record.annotations['strand2'] is '+': continue
		if record.annotations['exonIndex1'].startswith('0') and record.annotations['strand1'] is '-': continue
		
		# Filter out acceptor genes that include 3'-most exon for - strand acceptors (START codon),
		# or donor genes that end (start) with 3'-most exon (STOP codon) for + strand.
		(transA, transB) = record.id.split('-')
		if record.annotations['strand1'] == "+":
			record.annotations['geneAjunc'] = max(record.annotations['exonEnds1'].split(','))
			if checkAgainstRefFlat(transA,record.annotations['geneAjunc'] ,start=False) is None: continue
		else:
			record.annotations['geneAjunc'] = max(record.annotations['exonStarts1'].split(','))
			if checkAgainstRefFlat(transA,record.annotations['geneAjunc'] ,start=True) is None:  continue
		if record.annotations['strand2'] == "+":
			record.annotations['geneBjunc'] = min(record.annotations['exonStarts2'].split(','))
			if checkAgainstRefFlat(transB,record.annotations['geneBjunc'] ,start=True) is None: continue
		else: 
			record.annotations['geneBjunc'] = min(record.annotations['exonEnds2'].split(','))
			if checkAgainstRefFlat(transB,record.annotations['geneBjunc'] ,start=False) is None: continue

		
		# Transcript length filter
		if len(record) > minLength:
			filteredFusimEvents.append(record)

	if len(filteredFusimEvents) is not numEvents:
		print('Number of filtered events: '+len(filteredFusimEvents))
		
		
	with open(simName+'_filtered.gtf', 'w') as gtf, open(simName+'_filtered.bedpe', 'w') as bedpe:

		for record in filteredFusimEvents:

			# Write accepted events to GTF and BEDPE
			writeGTF(record_id=record.id,record_len=len(record),GTFfh=gtf)
			writeBEDPE(record.annotations,record.id,bedpe)
		
			# Write out separate fasta files per seq
			fastaFilename = '.'.join([record.id, 'fasta'])
			allFastaFilenames += fastaFilename+','
			numWritten = SeqIO.write(record, fastaFilename, "fasta")
	gtf.close()
	bedpe.close()

 	return(allFastaFilenames)
	
	

def makeFusionReference(fastaList, simName, numEvents, memory="2G", cores=4, disk="1M"):
	'''Runs RSEM to make reference for fusion events.'''
	
	cmd = ' '.join(['rsem-prepare-reference --gtf', simName+'_filtered.gtf', '--star --num-threads 4', fastaList.rstrip(','), '_'.join([simName, str(numEvents), 'ev'])])
	print(cmd)
	subprocess.call(cmd, shell=True)



###############
# Other functions
###############

def getHGNC(ENSG):
	r = requests.get("http://grch37.rest.ensembl.org/xrefs/id/"+ENSG+"?external_db=HGNC", headers={ "Content-Type" : "application/json"})
	if not r.ok:
		r.raise_for_status()
		return('HGNClookupFailed')
	decoded = r.json()[0]
	return(decoded['display_id'])


	
def parseFusimFastaDescription(record_description):
	""""Given a SeqRecord description string, returns the constituent items as a dictionary.

e.g. "fusionType=hybrid" -> annotations['fusionType']='hybrid'
"""
	annotations = dict()
	for item in record_description.strip().split():
		if item.startswith('ref'): continue
		k,v = item.split('=')
		annotations[k] = v
	return(annotations)



def parseFusimTxtOutput(fusimFilename):
	""""

e.g. 
"""
	additionalAnnotations = dict()
	with open(fusimFilename, 'r') as fusim:
		geneA = None # Parsing two lines at a time...
		for line in fusim:
			if not line.startswith('ENS'): continue
			if geneA is None:
				geneA = line.strip().split()
			else:
				geneB = line.strip().split()
				assert geneA[0] == geneB[0]
				tName = '-'.join([geneA[2],geneB[2]])
				additionalAnnotations[tName] = {'exonBases1':geneA[6],'exonStarts1':geneA[8],'exonEnds1':geneA[9],'exonBases2':geneB[6],'exonStarts2':geneB[8],'exonEnds2':geneB[9]}
				geneA = None
	fusim.close()
	return(additionalAnnotations)




def writeGTF(record_id,record_len,GTFfh):
	""""Writes a GTF entry for a gene fusion.

	"""
	seqname = record_id
	field9 = ' '.join(['gene_id "'+seqname+'";', 'gene_biotype "fusion";'])
	geneGtfLine = '\t'.join([seqname, 'fusim', 'gene', '1', str(record_len), '.', '+', '.', field9])
	GTFfh.write(geneGtfLine+'\n')
	field9 = ' '.join(['gene_id "'+seqname+'";', 'transcript_id "'+seqname+'";', 'gene_biotype "fusion";'])
	transcriptGtfLine = '\t'.join([seqname, 'fusim', 'transcript', '1', str(record_len), '.', '+', '.', field9])
	GTFfh.write(transcriptGtfLine+'\n')
	field9 = ' '.join(['gene_id "'+seqname+'";', 'transcript_id "'+seqname+'";', 'exon_number "1";', 'gene_biotype "fusion";'])
	gtfLine = '\t'.join([seqname, 'fusim', 'exon', '1', str(record_len), '.', '+', '.', field9])
	GTFfh.write(gtfLine+'\n')
	

	
def writeBEDPE(record_annotations,record_id,BEDPEfh):
	""""Writes a BEDPE entry for a gene fusion junction.
Note that FUSIM txt output is 1-based, although input refflat is 0-based.
	"""
	
	(transA, transB) = record_id.split('-')
	# subtract 1 from upstream position to be consistent with BED 0-based numbering	
	upstreamA = int(record_annotations['geneAjunc']) - 1
	upstreamB = int(record_annotations['geneBjunc']) - 1
#	bedpeTxt = '\t'.join([geneA[3], str(upstreamA), geneApos, geneB[3], str(upstreamB), geneBpos, geneA[0], '0', geneA[4], geneB[4], getHGNC(nameA), getHGNC(nameB)])
	bedpeTxt = '\t'.join([record_annotations['chrom1'], str(upstreamA), record_annotations['geneAjunc'], record_annotations['chrom2'], str(upstreamB), record_annotations['geneBjunc'], record_annotations['fusionGene'], '0', record_annotations['strand1'], record_annotations['strand2']])
	BEDPEfh.write(bedpeTxt+'\n')
	


def checkAgainstRefFlat(inENST,inPos,start=True):
	'''Check to see whether junction endpoints (1-based) match the refflat file (0-based). Also filter out certain cases of 
	unacceptable donor/acceptor candidates.
	'''
	refFlatBytes = subprocess.check_output(' '.join(['grep', inENST, referenceGeneModels]), shell = True)
	refflatLine = refFlatBytes.decode("utf-8").strip().split()
	if start:
		# subtract 1 to compare with 0-based refflat start positions
		return(re.search(str(int(inPos)-1), refflatLine[9]))
	else:
		if inPos == refflatLine[7]: return(None) # This line to remove cases where the stop codon exon is a donor 
													# for the plus strand, or the start codon is an acceptor for 
													# the minus strand. Checks for whether the position of interest
													# matches the end of the CDS.
		else: return(re.search(inPos, refflatLine[10]))



if __name__=="__main__":

	parser = argparse.ArgumentParser("Runs workflow to generate fusion events and truth file.")
#	Job.Runner.addToilOptions(parser)
	parser.add_argument('--numEvents', default=10, help='Number of filtered fusion events to generate.', type=int, required=False)
	parser.add_argument('--minLength', default=400, help='Minimum length of fusion transcript.', type=int, required=False)
	parser.add_argument("--simName", help="Prefix for the simulation filenames.", default='testSimulation', required=False)
	args = parser.parse_args()



# 	with open(os.environ['MODULESHOME']+'/init/python.py') as f:
# 		code = compile(f.read(), os.environ['MODULESHOME']+'/init/python.py', 'exec')
# 		exec(code)

#	exec(open(os.environ['MODULESHOME']+'/init/python.py').read())
#	module('load','rsem/1.2.8')

# 
# 	## Wrap jobs
# 	j1 = Job.wrapFn(runFusim, numEvents=args.numEvents, simName=args.simName)
# 	j2 = Job.wrapFn(filterFusionEvents, numEvents=args.numEvents, minLength=args.minLength, simName=args.simName)
# #	j3 = Job.wrapFn(makeFusionReference, fastaList=j2.rv(), simName=args.simName, numEvents=args.numEvents)
# 	
# 	## Specify order
# 	j1.addFollowOn(j2)
# #	j2.addFollowOn(j3)
# 	
# 	
# 	Job.Runner.startToil(j1, args)

	j1 = runFusim(numEvents=args.numEvents, simName=args.simName)
	j2 = filterFusionEvents(numEvents=args.numEvents, minLength=args.minLength, simName=args.simName)
