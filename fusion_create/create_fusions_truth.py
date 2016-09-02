#! /usr/bin/env python

from toil.job import Job
import subprocess
import argparse
import requests
import os

    

referenceGenome = '/external-data/Genome/genomes/Hsapiens_Ensembl_GRCh37/primary_nomask/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
referenceGeneModels = '/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/Hsapiens_Ensembl_v75.refflat.txt'

########################
## Workflow functions
########################

def runFusim(numEvents, simName, memory="2G", cores=1, disk="1M"):
	'''Runs Fusim to generate fusion events.'''
	
	fusimJarPath = '/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Scripts/fusim-0.2.2/fusim.jar'
	
	numEventsToSimulate = int(numEvents)*5
	cmd = ''.join(['java -jar ', fusimJarPath, ' --gene-model=', referenceGeneModels, ' --fusions=', str(numEventsToSimulate), ' --reference=', referenceGenome, ' --fasta-output=', simName+'.fasta', ' --text-output=', simName+'fusions.txt', ' --auto-correct-orientation --cds-only --keep-exon-boundary'])
	print cmd
	subprocess.call(cmd, shell = True)


def filterFusionEvents(numEvents, minLength, simName, memory="100M", cores=1, disk="1M"):
	'''Reads in Fusim output files and writes out a GTF, FASTAs, and BEDPE of events meeting min length criteria.'''
	
	filteredFusimEvents = dict()
	allFastaFilenames = ''


	# Go through the fasta file and filter it
	numberMeetingCriteria = 0
	inputFusim = simName+'.fasta'
	with open(inputFusim, 'r') as fusim:
		outfile = ''.join([simName, '_filtered', '.gtf'])
		with open(outfile, 'w') as gtf:
			for line in fusim:
				if numberMeetingCriteria >= numEvents: break
				elif line.startswith('>'):
					fullName = line
				else:
					length = len(line.strip())
					if length > minLength:
#						print fullName,line
						numberMeetingCriteria += 1
					
						# write out GTF
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
						
						genesFused = fullName.strip().split()[1].lstrip('fusionGene=')
						filteredFusimEvents[genesFused] = seqname
#						print ' '.join(['Adding to filteredFusimEvents dict:', genesFused, seqname])
						
						# write out fasta
						fastaFilename = '.'.join([seqname, 'fasta'])
						allFastaFilenames += fastaFilename+','
						fasta = open(fastaFilename,'w')
						fasta.write('>'+seqname+'\n')
						fasta.write(line)
						fasta.close()
		gtf.close()
	fusim.close()
						
	# Make truth file
	with open(simName+'fusions.txt', 'r') as fusim2:
	
		with open(simName+'_filtered.bedpe', 'w') as bedpe:

			geneA = None
			for line in fusim2:
		#		print line
				if not line.startswith('ENS'): continue
				if geneA is None:
					geneA = line.strip().split()
					if geneA[4] == "+":
						geneApos = max(geneA[9].split(','))
					else:
						geneApos = max(geneA[8].split(','))
		#				print '%s' % geneApos
				else:
					geneB = line.strip().split()
			
					# skip entries in this file that are not part of the filtered set
					if geneA[0] in filteredFusimEvents: 
						transcriptsJoint = '-'.join([geneA[2], geneB[2]])
						if filteredFusimEvents[geneA[0]] == transcriptsJoint:
							(nameA, nameB) = geneA[0].split('-')
							if geneB[4] == "+":
								geneBpos = min(geneB[8].split(','))
							else: 
								geneBpos = min(geneB[9].split(','))
			
							# subtract 1 from upstream position to be consistent with BED 0-based numbering	
							upstreamA = int(geneApos) - 1
							upstreamB = int(geneBpos) - 1
#							bedpeTxt = '\t'.join([geneA[3], str(upstreamA), geneApos, geneB[3], str(upstreamB), geneBpos, geneA[0], '0', geneA[4], geneB[4], getHGNC(nameA), getHGNC(nameB)])
							bedpeTxt = '\t'.join([geneA[3], str(upstreamA), geneApos, geneB[3], str(upstreamB), geneBpos, geneA[0], '0', geneA[4], geneB[4]])
							bedpe.write(bedpeTxt+'\n')
#							print '%s' % bedpeTxt
					geneA = None
					geneB = None
	
	return(allFastaFilenames)
	
	

def makeFusionReference(fastaList, simName, numEvents, memory="2G", cores=4, disk="1M"):
	'''Runs RSEM to make reference for fusion events.'''
	
	cmd = ' '.join(['rsem-prepare-reference --gtf', simName+'_filtered.gtf', '--star --num-threads 4', fastaList.rstrip(','), '_'.join([simName, str(numEvents), 'ev'])])
	print cmd
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



if __name__=="__main__":

	parser = argparse.ArgumentParser("Runs workflow to generate fusion events and truth file.")
	Job.Runner.addToilOptions(parser)
	parser.add_argument('--numEvents', default=10, help='Number of filtered fusion events to generate.', type=int, required=False)
	parser.add_argument('--minLength', default=400, help='Minimum length of fusion transcript.', type=int, required=False)
	parser.add_argument("--simName", help="Prefix for the simulation filenames.", default='testSimulation', required=False)
	args = parser.parse_args()

	execfile(os.environ['MODULESHOME']+'/init/python.py')
	module('load','rsem/1.2.8')


	## Wrap jobs
	j1 = Job.wrapFn(runFusim, numEvents=args.numEvents, simName=args.simName)
	j2 = Job.wrapFn(filterFusionEvents, numEvents=args.numEvents, minLength=args.minLength, simName=args.simName)
	j3 = Job.wrapFn(makeFusionReference, fastaList=j2.rv(), simName=args.simName, numEvents=args.numEvents)
	
	## Specify order
	j1.addFollowOn(j2)
	j2.addFollowOn(j3)
	
	
	Job.Runner.startToil(j1, args)
