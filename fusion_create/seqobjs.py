
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def readGenome(fasta):
   genome_dict = SeqIO.index(fasta, "fasta")
   print(len(genome_dict))
   return(genome_dict)
   
   
def makeFusionSeqObj(donorExonSeq,acceptorExonSeq,genomeObj):

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
   fusionObj.annotations['dJunction'] = 99
   fusionObj.annotations['aJunction'] = 99
   return(fusionObj)
   
   

def concatExonSeq(exonList,genomeObj):

   chrom = exonList[0].qualifiers['seqid'][0]
   fusion_seq = ''
   if exonList[0].strand == 1:
      for exon in exonList:
         fusion_seq = fusion_seq + exon.extract(genomeObj[chrom])
         printExonSF(exon)
   else:
      exonList.reverse()
      for exon in exonList:
         fusion_seq = fusion_seq + exon.extract(genomeObj[chrom])
         printExonSF(exon)


def printExonSF(exonSF):

   print(exonSF.qualifiers['exon_id'], exonSF.location.start, exonSF.location.end, exonSF.qualifiers['exon_number'], exonSF.id)




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
   # subtract 1 from upstream position to be consistent with BED 0-based numbering   
   upstreamA = int(record.annotations['dJunction']) - 1
   upstreamB = int(record.annotations['aJunction']) - 1
   bedpeTxt = '\t'.join([record.annotations['dChrom'], str(upstreamA), record.annotations['dJunction'], record.annotations['aChrom'], str(upstreamB), record.annotations['aJunction'], record.annotations['fusionGene'], '0', record.annotations['dStrand'], record.annotations['aStrand']])
   BEDPEfh.write(bedpeTxt+'\n')