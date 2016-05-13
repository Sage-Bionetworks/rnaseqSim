#! /usr/bin/env python

import subprocess
import argparse
import os
import random

    


########################
## Workflow functions
########################

def generateReads(model, isoV, simName, fusRef, fusV, simReads, otherReads, memory="2G", cores=1, disk="15G"):
	'''Runs Fusim to generate fusion events.'''


	cmd = ' '.join(['rsem-simulate-reads /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/diploid_reference_genome/STAR/GRCh37v75_diploid', model, isoV, '0.066', str(otherReads), simName+'_diploid'])
	print cmd
#	subprocess.call(cmd, shell = True)

	cmd = ' '.join(['rsem-simulate-reads', fusRef, model, fusV, '0.066', str(int(simReads)), simName+'_fusions'])
	print cmd
#	subprocess.call(cmd, shell = True)



def postProcessReads(simName, totalReads, simReads, memory="100M", cores=1, disk="15G"):
	'''Changes read names and merges read files.'''
	
	insertPosition = random.sample(xrange(int(totalReads)),int(simReads))
	insertPosition.sort(reverse=True)
	print '%s' % insertPosition
	insertPosition2 = list(insertPosition)
	
	diploidR1 = simName+'_diploid_1.fq'
	diploidR2 = simName+'_diploid_2.fq'
	fusionR1 = simName+'_fusions_1.fq'
	fusionR2 = simName+'_fusions_2.fq'
	
	renameAndMerge(diploid = diploidR1, fusion = fusionR1, inserts = insertPosition, keyFileName = simName+'_readKey_1.txt', fastqFileName = simName+'_merged_1.fq')
	
	renameAndMerge(diploid = diploidR2, fusion = fusionR2, inserts = insertPosition2, keyFileName = simName+'_readKey_2.txt', fastqFileName = simName+'_merged_2.fq')
	



###############
# Other functions
###############

def renameAndMerge(diploid, fusion, inserts, keyFileName, fastqFileName):

	keyFH = open(keyFileName, 'w')
	fqFH = open(fastqFileName, 'w')
	
	fusionFH = open(fusion, 'r')

	with open(diploid, 'r') as diploidFQ:
		lineCounter = 0
		fusRead = ''
		diploidRead = ''
		insertsCurrent = inserts.pop()
		for line in diploidFQ:
			# Check for complete FASTQ record
			if lineCounter % 4 == 0:
				# Check for whether to insert one of the fusion reads
				print lineCounter, lineCounter/4, insertsCurrent
				if lineCounter/4 == insertsCurrent:
					for i in xrange(4):
						tmp = ''.join([fusRead,fusionFH.readline()])
						fusRead = tmp
					ifrReturns = insertFusionRead(fusionRead=fusRead, inserts=inserts, lC=lineCounter, dR=diploidRead, keyFH=keyFH, fqFH=fqFH, fusionFH=fusionFH)
					print ifrReturns
					insertsCurrent = ifrReturns[0]
					lineCounter = ifrReturns[1]
					fusRead = ''
				elif diploidRead == '':
					pass
				else: 
					writeFastQ(readString=diploidRead, readNumber=lineCounter/4, keyFileHandle=keyFH, fastqFileHandle=fqFH)
				diploidRead = line
			else:
				tmp = ''.join([diploidRead, line])
				diploidRead = tmp
			lineCounter += 1
		writeFastQ(readString=diploidRead, readNumber=lineCounter/4, keyFileHandle=keyFH, fastqFileHandle=fqFH)
	
	keyFH.close()
	fqFH.close()
	

def insertFusionRead(fusionRead, inserts, lC, dR, keyFH, fqFH, fusionFH):

	writeFastQ(readString=fusionRead, readNumber=lC/4, keyFileHandle=keyFH, fastqFileHandle=fqFH)
	if len(inserts) > 0:
		insertsCr = inserts.pop()
	else:
		insertsCr = -1
	lC += 4
	# check that the next read is not also a fusion read
	if not insertsCr == lC/4:
		writeFastQ(readString=dR, readNumber=lC/4, keyFileHandle=keyFH, fastqFileHandle=fqFH)
		return([insertsCr, lC])
	else:
		fusRead = ''
		for i in xrange(4):
			tmp = ''.join([fusRead, fusionFH.readline()])
			fusRead = tmp
		insertFusionRead(fusionRead=fusRead, inserts=inserts, lC=lC, dR=dR, keyFH=keyFH, fqFH=fqFH, fusionFH=fusionFH)


def writeFastQ(readString, readNumber, keyFileHandle, fastqFileHandle):
	
	lines = readString.split('\n')
	keyFileHandle.write('\t'.join([str(readNumber), lines[0], '\n']))
	fastqFileHandle.write('\n'.join(['@'+str(readNumber),lines[1],lines[2],lines[3]+'\n']))




if __name__=="__main__":

	parser = argparse.ArgumentParser("Runs workflow to generate fusion reads files.")
	parser.add_argument('--totalReads', default=5e6, help='Total number of reads to generate.', type=int, required=False)
	parser.add_argument('--numSimReads', default=5e5, help='Total number of simulated reads to generate.', type=int, required=False)
	parser.add_argument("--simName", help="Prefix for the simulation filenames.", default='testSimulation', required=False)
	parser.add_argument("--RSEMmodel", help="Model file from RSEM alignment.", required=True)
	parser.add_argument("--isoformTPM", help="File of isoform TPM values to simulate.", required=True)
	parser.add_argument("--fusionTPM", help="File of fusion TPM values to simulate.", required=True)
	parser.add_argument("--fusRef", help="Path to fusion RSEM-format reference.", required=True)
	args = parser.parse_args()

	execfile(os.environ['MODULESHOME']+'/init/python.py')
	module('load','rsem/1.2.8')
	module('load','star/2.4.2a')



	## Wrap jobs
	generateReads(model=args.RSEMmodel, isoV=args.isoformTPM, simName=args.simName, fusRef=args.fusRef, fusV=args.fusionTPM, simReads=args.numSimReads, otherReads=args.totalReads-args.numSimReads)
	postProcessReads(simName=args.simName, totalReads=args.totalReads, simReads=args.numSimReads)
	
