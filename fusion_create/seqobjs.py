
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
   fusionObj = SeqRecord(Seq(str(donor_seq)+str(acceptor_seq)), id=fusionId)
   fusionObj.annotations['dStrand'] = donorExonSeq[0].strand
   fusionObj.annotations['aStrand'] = acceptorExonSeq[0].strand
   fusionObj.annotations['dChrom'] = donorChrom
   fusionObj.annotations['aChrom'] = acceptorChrom
   fusionObj.annotations['dJunction'] = dJunc
   fusionObj.annotations['aJunction'] = aJunc
   return(fusionObj)
   
   

def concatExonSeq(exonList,genomeObj):

   chrom = exonList[0].qualifiers['seqid'][0]
   fusion_seq = ''
   if exonList[0].strand == '-1':
      exonList.reverse()
   for exon in exonList:
      tmp_seq = exon.extract(genomeObj[chrom])
      printExonSF(exon,tmp_seq)
      fusion_seq = fusion_seq + tmp_seq
   return(fusion_seq)


def printExonSF(exonSF,seq):

   print(exonSF.qualifiers['exon_id'], exonSF.location.start, exonSF.location.end, exonSF.qualifiers['exon_number'], exonSF.id,str(seq))




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