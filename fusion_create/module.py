#!/usr/bin/env python

import subprocess
import argparse
import os
#import gzip
import gffutils
import fusions
import seqobjs
import random
from Bio import SeqIO

def run_module(genome_file, 
               gtf_file, 
               numEvents, 
               simName, 
               mid_exon_prob,
               mid_exon_min_size,
               mid_exon_min_cleaved,
               add_sense_antisense_fusions = False,
               add_exon_duplications_and_deletions = False,
               add_fusion_events_in_UTR = False):
    
    # Converting GTF file into a database
    database_filename = '.'.join([os.path.basename(gtf_file).rstrip('.gtf'), 'sqlite3'])     
    dbPath = os.path.join(os.path.dirname(gtf_file),database_filename)
    print(dbPath)

    if os.path.isfile(dbPath):
        # Connect to an already-existing db
        db = gffutils.FeatureDB(dbPath)
    else:
        # Or, create a new one
#        db = gffutils.create_db(gtf_file,database_filename)
        db = gffutils.create_db(gtf_file,database_filename,disable_infer_genes=True, disable_infer_transcripts=True)

    # Filter the gene types to consider, e.g. protein-coding
    protein_coding_genes = list()
    allGenesIter = db.features_of_type("gene")
    for item in allGenesIter:
        if item['gene_biotype'][0] == 'protein_coding':
            protein_coding_genes.append(item.id)

    # Get the number of genes available after filtering     
    print(' '.join(['Number of protein-coding genes:', str(len(protein_coding_genes))]))
    
    hg19 = seqobjs.readGenome(genome_file)
    fastaFilenames = list()    

    with open(''.join([simName, '.gtf']),'w') as gtf, open(''.join([simName, '_filtered.bedpe']),'w') as bedpe:
    # Get fusion events as tuples of Bio.Seq objects
    # TODO: Simplify return objects, possibly returning one event at a time instead of list
        for fusion_event in fusions.getRandomFusions(
            db = db, names = protein_coding_genes, num = numEvents):
            fusion_event.create_breakages(
                mid_exon_prob, mid_exon_min_size, mid_exon_min_cleaved)
            fObj = seqobjs.makeFusionSeqObj(
                donorExonSeq = fusion_event.get_donor_exons(), 
                acceptorExonSeq = fusion_event.get_acceptor_exons(),
                dJunc = fusion_event.get_donor_junction(),
                aJunc = fusion_event.get_acceptor_junction(),
                genomeObj = hg19)
            print(len(fObj))
            seqobjs.writeGTF(fObj,gtf)
            seqobjs.writeBEDPE(fObj,bedpe)
            SeqIO.write(fObj, ''.join([fObj.id,'.fasta']), "fasta")
            fastaFilenames.append(''.join([fObj.id,'.fasta']))
    
    return(fastaFilenames)
    
    

def makeFusionReference(fastaList, simName, numEvents):
   '''Runs RSEM to make reference for fusion events.'''
   
   cmd = ' '.join(['rsem-prepare-reference --gtf', simName+'.gtf', '--star --num-threads 4', ','.join(fastaList), '_'.join([simName, str(numEvents), 'ev'])])
   print(cmd)
   subprocess.call(cmd, shell=True)


    
if __name__ == '__main__':


    parser = argparse.ArgumentParser("Runs workflow to generate fusion events and truth file.")
    parser.add_argument('--genome', default='/external-data/Genome/genomes/Hsapiens_Ensembl_GRCh37/primary_nomask/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa', help='Reference Genome.', type=str, required=False)
    parser.add_argument('--gtf', default='/external-data/Genome/gene_models/Hsapiens_Ensembl_v75_refonly.gtf', help='Gene models in GTF format.', type=str, required=False)
    parser.add_argument('--numEvents', default=5, help='Number of filtered fusion events to generate.', type=int, required=False)
    parser.add_argument('--minLength', default=400, help='Minimum length of fusion transcript.', type=int, required=False)
    parser.add_argument("--simName", help="Prefix for the simulation filenames.", default='testSimulation', required=False)
    parser.add_argument(
        "--seed", 
        help = "Seed number to use for RSEM read simulation.", 
        type = int, 
        required = False, 
        default = None)
    # mid exon fusion parameters ----------------------------------------------                  
    parser.add_argument(
        '--mid_exon_prob', 
        default = 0.0, 
        help = 'Probability of mid exon breakage per donor/acceptor',
        type = float,
        required = False)
    parser.add_argument(
        '--mid_exon_min_size', 
        default = 1, 
        help = 'Min exon size after mid-exon breakage',
        type = int,
        required = False)
    parser.add_argument(
        '--mid_exon_min_cleaved', 
        default = 1, 
        help = 'Min bases removed per mid-exon breakage',
        type = int,
        required = False)

    
    args = parser.parse_args()
    
    # set seed to seed arument
    if isinstance(args.seed, (int, long)):
        random.seed(args.seed)
                     
    fastaFN = run_module(args.genome, 
                         args.gtf,
                         args.numEvents, 
                         args.simName,
                         args.mid_exon_prob,
                         args.mid_exon_min_size,
                         args.mid_exon_min_cleaved)
                         
    makeFusionReference(fastaList = fastaFN, 
                        simName = args.simName, 
                        numEvents = args.numEvents)
