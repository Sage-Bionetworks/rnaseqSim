
import sys
import gzip
import gtf
import random
def run_module(genome_file, gtf_file):
    in_handle = gzip.GzipFile(gtf_file)
    g = gtf.GTF()
    g.read(in_handle)

    transcripts = g.transcript_collect()

    for i in range(100):
        fgene = [ 
            transcripts[random.choice(transcripts.keys())],
            transcripts[random.choice(transcripts.keys())]
        ]
        
        for j in fgene:
            exons = sorted( list(k for k in j if 'exon_number' in k.attribute), key=lambda x: int(x.attribute['exon_number']) )
            
            print j[0].attribute['gene_name'], ",".join( 
                ("%s(%s:%s-%s:%s)" % (k.attribute['exon_number'], k.reference, k.start, k.end, k.strand) for k in exons)
                )

        print "-=-=-"
    


if __name__ == '__main__':
    run_module(sys.argv[1], sys.argv[2])