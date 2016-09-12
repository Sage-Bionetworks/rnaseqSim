
import sys
import os
#import gzip
import gffutils


def run_module(genome_file, gtf_file):

   # Converting GTF file into a database
   database_filename = '.'.join([os.path.basename(gtf_file).rstrip('.gtf'), 'sqlite3'])    

   if os.path.isfile(database_filename):
      # Connect to an already-existing db
      db = gffutils.FeatureDB(database_filename)
   else:
      # Or, create a new one
      db = gffutils.create_db(gtf_file,database_filename)


   # Filter the gene types to consider, e.g. protein-coding
   protein_coding_genes = list()
   allGenesIter = db.features_of_type("gene")
   for item in allGenesIter:
      if item['gene_biotype'][0] == 'protein_coding':
         protein_coding_genes.append(item.id)

   # Get the number of genes available after filtering    
   print(' '.join(['Number of protein-coding genes:', str(len(protein_coding_genes))]))

   # Get the transcripts for a gene, returns generator
   ENSG = 'ENSG00000239713'
   transcripts = db.children(ENSG,featuretype = "transcript")
   # Put the protein-coding transcripts in a list
   protein_coding_transcripts = list()
   for item in transcripts:
      if item.source == 'protein_coding':
         protein_coding_transcripts.append(item)
   print(' '.join(['Number of protein-coding transcripts for gene', ENSG+':', str(len(protein_coding_transcripts))]))

if __name__ == '__main__':
   run_module(sys.argv[1], sys.argv[2])