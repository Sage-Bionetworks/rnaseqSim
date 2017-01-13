#! /usr/bin/env python

import subprocess
import argparse
import os
import random
import shutil
import synapseclient

    


########################
## Workflow functions
########################

def generateReads(modelID, isoV, simName, fusRef, fusV, simReads, dipGenome, otherReads, memory="2G", cores=1, disk="15G"):
	'''Runs Fusim to generate fusion events.'''


	model = syn.get(modelID, downloadLocation = os.getcwd())
	cmd = ' '.join(['rsem-simulate-reads', dipGenome, model.path, isoV, '0.066', str(otherReads), simName+'_diploid'])
	print cmd
	subprocess.call(cmd, shell = True)

	cmd = ' '.join(['rsem-simulate-reads', fusRef, model.path, fusV, '0.066', str(int(simReads)), simName+'_fusions'])
	print cmd
	subprocess.call(cmd, shell = True)



def postProcessReads(simName, totalReads, simReads, memory="100M", cores=1, disk="15G"):
	'''Changes read names and merges read files.'''
	
	insertPosition = random.sample(xrange(int(totalReads)),int(simReads))
	insertPosition.sort(reverse=False)
#	print '%s' % insertPosition
	
	diploidR1 = simName+'_diploid_1.fq'
	diploidR2 = simName+'_diploid_2.fq'
	fusionR1 = simName+'_fusions_1.fq'
	fusionR2 = simName+'_fusions_2.fq'
	
	os.mkdir('tmp')
	
	renameAndMerge(diploid = diploidR1, fusion = fusionR1, inserts = insertPosition, keyFileName = simName+'_readKey_1.txt', fastqFileName = simName+'_merged_1.fq')
	cmd = ' '.join(['sort -T tmp -k1,1', simName+'_merged_1.fq', '| sed -e \'s/\\t/\\n/g\' - | gzip >', simName+'_mergeSort_1.fq.gz'])
	print cmd
	subprocess.call(cmd, shell = True)
	
	
	renameAndMerge(diploid = diploidR2, fusion = fusionR2, inserts = insertPosition, keyFileName = simName+'_readKey_2.txt', fastqFileName = simName+'_merged_2.fq')
	cmd = ' '.join(['sort -T tmp -k1,1', simName+'_merged_2.fq', '| sed -e \'s/\\t/\\n/g\' - | gzip >', simName+'_mergeSort_2.fq.gz'])
	print cmd
	subprocess.call(cmd, shell = True)
	
	shutil.rmtree('tmp')
	
	

def makeIsoformsTruth(simName, memory="100M", cores=1, disk="50M"):
	'''Takes in RSEM isoforms.results file from diploid simulation.'''
	
	cmd = ' '.join(['sort', simName+'_diploid.sim.isoforms.results >', simName+'_diploid.sim.isoforms.results_sorted'])
	print cmd
	subprocess.call(cmd, shell = True)
	
	truthFH = open(simName+'_isoforms_truth.txt', 'w')
	with open(simName+'_diploid.sim.isoforms.results_sorted', 'r') as rsem:
		line1 = None
		for line in rsem:
			if not line.startswith('ENST'): continue
			if line1 is None:
				line1 = line
			else:
				# check that transcript ids match
				line1v = line1.strip().split()
				line2v = line.strip().split()
				if not line1v[0].split('-')[0] == line2v[0].split('-')[0]:
					print 'Line matching off %s %s' % (line1v[0], line2v[0])
					break
				else:
					summedTPM = float(line1v[5]) + float(line2v[5])
					truthFH.write('%s\t%f\n' % (line1v[0].split('-')[0], summedTPM))
				line1 = None
	truthFH.close()



###############
# Other functions
###############

def renameAndMerge(diploid, fusion, inserts, keyFileName, fastqFileName):

	keyFH = open(keyFileName, 'w')
	fqFH = open(fastqFileName, 'w')
	
	# Write out fusion reads with new numbers
	iC = 0 # inserts counter
	read = ''
	with open(fusion, 'r') as fusionFH:
		lC = 0
		for line in fusionFH:
			if lC == 0:
				header = '\t'.join([str(inserts[iC]),line.strip()])	
				read = '@'+str(inserts[iC])
				iC += 1			
			elif lC % 4 == 0 and lC > 0:
				writeFastQ(header, read, keyFileHandle=keyFH, fastqFileHandle=fqFH) 				
				read = '@'+str(inserts[iC])
				header = '\t'.join([str(inserts[iC]),line.strip()])	
				iC += 1		
			else:
				tmp = '\t'.join([read, line.strip()])
				read = tmp	
			lC +=1
		writeFastQ(header, read, keyFileHandle=keyFH, fastqFileHandle=fqFH) 				

	# Write out regular reads with new numbers
	rC = 0 # read counter
	read = ''
	lC = 0 # line counter
	iC = 0 # inserts counter
	with open(diploid, 'r') as diploidFH:
		for line in diploidFH:
			# skip read numbers that were already assigned to fusion reads
			if iC < len(inserts):
				while rC == inserts[iC]:
					rC += 1
					iC += 1
					if iC >= len(inserts): break
			# print reads with new numbers
			if lC % 4 == 0:		
				if lC > 0:
					writeFastQ(header, read, keyFileHandle=keyFH, fastqFileHandle=fqFH) 								
				header = '\t'.join([str(rC),line.strip()])	
				read = '@'+str(rC)
				rC += 1			
			else:
				tmp = '\t'.join([read, line.strip()])
				read = tmp	
			lC +=1
		# print the last read
		writeFastQ(header, read, keyFileHandle=keyFH, fastqFileHandle=fqFH) 				

	
	keyFH.close()
	fqFH.close()
	

def writeFastQ(header, readString, keyFileHandle, fastqFileHandle):
	
	keyFileHandle.write(header+'\n')
	fastqFileHandle.write(readString+'\n')




if __name__=="__main__":

	parser = argparse.ArgumentParser("Runs workflow to generate fusion reads files.")
	parser.add_argument('--totalReads', default=5e6, help='Total number of reads to generate.', type=int, required=False)
	parser.add_argument('--numSimReads', default=5e5, help='Total number of simulated reads to generate.', type=int, required=False)
	parser.add_argument("--simName", help="Prefix for the simulation filenames.", default='testSimulation', required=False)
	parser.add_argument("--RSEMmodel", help="Model file from RSEM alignment.", required=True)
	parser.add_argument("--isoformTPM", help="File of isoform TPM values to simulate.", required=True)
	parser.add_argument("--fusionTPM", help="File of fusion TPM values to simulate.", required=True)
	parser.add_argument("--fusRef", help="Path to fusion RSEM-format reference.", required=True)
        parser.add _argument("--dipGenome", help="File of the diploid genome.", required=True)
	args = parser.parse_args()
	
	syn = synapseclient.login()

	## Wrap jobs
	generateReads(modelID=args.RSEMmodel, isoV=args.isoformTPM, simName=args.simName, fusRef=args.fusRef, fusV=args.fusionTPM, simReads=args.numSimReads, dipGenome=args.dipGenome, otherReads=args.totalReads-args.numSimReads)
	postProcessReads(simName=args.simName, totalReads=args.totalReads, simReads=args.numSimReads)
	makeIsoformsTruth(simName=args.simName)
	
