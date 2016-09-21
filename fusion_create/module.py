
import sys
import os
#import gzip
import gffutils
import fusions
import seqobjs


def run_module(genome_file, gtf_file):

   # Converting GTF file into a database
   database_filename = '.'.join([os.path.basename(gtf_file).rstrip('.gtf'), 'sqlite3'])    

   if os.path.isfile(database_filename):
      # Connect to an already-existing db
      db = gffutils.FeatureDB(database_filename)
   else:
      # Or, create a new one
#      db = gffutils.create_db(gtf_file,database_filename)
      db = gffutils.create_db(gtf_file,database_filename,disable_infer_genes=True, disable_infer_transcripts=True)

   # Filter the gene types to consider, e.g. protein-coding
   protein_coding_genes = list()
   allGenesIter = db.features_of_type("gene")
   for item in allGenesIter:
      if item['gene_biotype'][0] == 'protein_coding':
         protein_coding_genes.append(item.id)

   # Get the number of genes available after filtering    
   print(' '.join(['Number of protein-coding genes:', str(len(protein_coding_genes))]))
   
   hg19 = seqobjs.readGenome(sys.argv[1])   
   for item in hg19.keys():
      print(item)

   # Get fusion events as tuples of Bio.Seq objects
   fusionEvents = fusions.getRandomFusions(db=db, names=protein_coding_genes)
   for event in fusionEvents:
      donorSeq,acceptorSeq = fusions.convertToSeqObj(event)
#      print('donor',donorSeq.qualifiers['seqid'][0], donorSeq.qualifiers['junctionExonNum'][0],donorSeq.qualifiers['source'][0],donorSeq.strand)
#      print('acceptor',acceptorSeq.qualifiers['seqid'][0], acceptorSeq.qualifiers['junctionExonNum'][0],acceptorSeq.qualifiers['source'][0],acceptorSeq.strand)
      fObj = seqobjs.makeFusionSeqObj(donorSeq,acceptorSeq,hg19)
      print(fObj.id)    
      

   
   
   
if __name__ == '__main__':
   run_module(genome_file=sys.argv[1], gtf_file=sys.argv[2])