
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def readGenome(fasta):
    genome_dict = SeqIO.index(fasta, "fasta")
    print(len(genome_dict))
    return(genome_dict)
    
    
def makeFusionSeqObj(donorExonSeq,acceptorExonSeq,dJunc,aJunc,genomeObj):

    dName = donorExonSeq[0].qualifiers['transcript_id'][0]
    aName = acceptorExonSeq[0].qualifiers['transcript_id'][0]
    fusionId = '-'.join([dName, aName])
    donorChrom = donorExonSeq[0].qualifiers['seqid'][0]
    acceptorChrom = acceptorExonSeq[0].qualifiers['seqid'][0]
    
    # Now make new Seq Record object containing concatenated fusion sequence    
    donor_seq = concatExonSeq(donorExonSeq,genomeObj)
    acceptor_seq = concatExonSeq(acceptorExonSeq,genomeObj)
    fusion_rec = SeqRecord(donor_seq+acceptor_seq, id=fusionId)
    fusion_rec.annotations['dStrand'] = donorExonSeq[0].strand
    fusion_rec.annotations['aStrand'] = acceptorExonSeq[0].strand
    fusion_rec.annotations['dChrom'] = donorChrom
    fusion_rec.annotations['aChrom'] = acceptorChrom
    fusion_rec.annotations['dJunction'] = dJunc
    fusion_rec.annotations['aJunction'] = aJunc
    return(fusion_rec)
    
    

def concatExonSeq(exonList,genomeObj):
    """"Returns Seq object containing concatenated seq of exons in the input list, thus, a partial transcript."""

    chrom = exonList[0].qualifiers['seqid'][0]
    fusion_seq = ''
    if exonList[0].strand == '-1':
        exonList.reverse()
    for exon in exonList:
        tmp_seq = exon.extract(genomeObj[chrom])
        printExon_seqfeat(exon,tmp_seq.seq)
        fusion_seq += tmp_seq
    return(fusion_seq.seq)


def printExon_seqfeat(exon_seqfeat,seq):

    print(exon_seqfeat.qualifiers['exon_id'], exon_seqfeat.location.start, exon_seqfeat.location.end, exon_seqfeat.qualifiers['exon_number'], exon_seqfeat.id)




def writeGTF(record,GTFfh):
    """"Writes a GTF entry for a gene fusion."""

    seqname = record.id
    field9 = ' '.join(['gene_id "'+seqname+'";', 'gene_biotype "fusion";'])
    geneGtfLine = '\t'.join([seqname, 'simulated', 'gene', '1', str(len(record)), '.', '+', '.', field9])
    GTFfh.write(geneGtfLine+'\n')
    field9 = ' '.join(['gene_id "'+seqname+'";', 'transcript_id "'+seqname+'";', 'gene_biotype "fusion";'])
    transcriptGtfLine = '\t'.join([seqname, 'simulated', 'transcript', '1', str(len(record)), '.', '+', '.', field9])
    GTFfh.write(transcriptGtfLine+'\n')
    field9 = ' '.join(['gene_id "'+seqname+'";', 'transcript_id "'+seqname+'";', 'exon_number "1";', 'gene_biotype "fusion";'])
    gtfLine = '\t'.join([seqname, 'simulated', 'exon', '1', str(len(record)), '.', '+', '.', field9])
    GTFfh.write(gtfLine+'\n')
    

    
def writeBEDPE(record,BEDPEfh):
    """"Writes a BEDPE entry for a gene fusion junction."""
    
    (transA, transB) = record.id.split('-')
    # add 1 to position to be consistent with BED 0-based numbering    
    downstreamA = int(record.annotations['dJunction']) + 1
    downstreamB = int(record.annotations['aJunction']) + 1
    bedpeTxt = '\t'.join([str(record.annotations['dChrom']), str(record.annotations['dJunction']), str(downstreamA), str(record.annotations['aChrom']), str(record.annotations['aJunction']), str(downstreamB), record.id, '0', str(record.annotations['dStrand']), str(record.annotations['aStrand'])])
    BEDPEfh.write(bedpeTxt+'\n')