
import sys
import os
#import gzip
import gffutils
import fusions


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
   
   # Get fusion events
   fusionEvents = fusions.getRandomFusions(db=db, names=protein_coding_genes, num=10)
   for event in fusionEvents:
      print(event)
   
   
if __name__ == '__main__':
   run_module(genome_file=sys.argv[1], gtf_file=sys.argv[2])