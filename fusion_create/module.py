
import sys
#import gzip
import gffutils


def run_module(genome_file, gtf_file):

	# Converting GTF file into a database
    database_filename = '.'.join([gtf_file.rstrip('gtf'), 'sqlite3'])    

    try:
	    # Connect to an already-existing db
	    db = gffutils.FeatureDB(database_filename)
	except FileNotFoundError:
		# Or, create a new one
		db = gffutils.create_db(gtf_file,database_filename)
    
    
    # Filter the gene types to consider, e.g. protein-coding
    protein_coding_genes = list()
    allGenesGen = db.somefunction(,)
	for item in allGenesGen:
		if item['gene_biotype'] is 'protein_coding':
			protein_coding_genes.append(item['Name'])

	# Get the number of genes available after filtering    
    number_genes = 

	# Get the transcripts for a gene, returns generator
	transcripts = db.children(ENSG,featuretype = "transcript")
	protein_coding_transcripts = list()
	for item in transcripts:
		

if __name__ == '__main__':
    run_module(sys.argv[1], sys.argv[2])