#! /bin/bash


# fusim prep
java -Xmx6G -jar fusim-0.2.2/fusim.jar --convert -i /external-data/Genome/gene_models/Hsapiens_gencode.v23.primary_assembly.annotation.gtf -o ../../Data/fusion-simdata-rnaf/Hsapiens_gencode.v23.primary_assembly.refflat.txt


java -jar fusim.jar --gene-model=/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/Hsapiens_gencode.v24.primary_assembly.refflat.txt --fusions=40 --reference=/external-data/Genome/genomes/Hsapiens_Gencode_GRCh38/primary/GRCh38.primary_assembly.genome.fa --fasta-output=fusions_2gene40exon.fasta --text-output=fusions_2gene40exon.txt --cds-only --keep-exon-boundary



module load star/2.4.2a
module load rsem

# this command doesn't work
rsem-prepare-reference --transcript-to-gene-map /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Results/fusions/rsem/gencode-fusion40-map-long.txt --star --star-sjdboverhang 74 --num-threads 4  /external-data/Genome/transcriptomes/Hsapiens_Gencode_v24/gencode.v24.transcripts.fa,/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Results/fusions/FUSIM/2genes_exonboundaries_40events/fusions_2gene40exon.fasta genv24_40event

# this command works
rsem-prepare-reference --gtf /external-data/Genome/gene_models/Hsapiens_Ensembl_v75.gtf --star --num-threads 4  /external-data/Genome/genomes/Hsapiens_Ensembl_GRCh37/primary_nomask/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa GRCh37v75_STAR




module load bowtie

~/Software/RSEM-1.2.28/rsem-prepare-reference --transcript-to-gene-map /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Results/fusions/rsem/gencode-fusion40-map-long.txt --bowtie2 --polyA  /external-data/Genome/transcriptomes/Hsapiens_Gencode_v24/gencode.v24.transcripts.fa,/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Results/fusions/FUSIM/2genes_exonboundaries_40events/fusions_2gene40exon.fasta genv24_40event


UrQt --in 110715_UNC14-SN744_0143_BC01UKABXX.8_1.fastq --inpair 110715_UNC14-SN744_0143_BC01UKABXX.8_2.fastq --out 110715_UNC14-SN744_0143_BC01UKABXX.8_1_trim.fastq --outpair 110715_UNC14-SN744_0143_BC01UKABXX.8_2_trim.fastq --gz --m 4 --min_read_size 40 --phred 33 --t 20 --buffer 500000 


rsem-calculate-expression --num-threads 4 --bowtie2  --estimate-rspd --paired-end /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/KIRC_data_and_RSEM_model/21a23647-5bce-4663-9c63-b1456c21503e/110715_UNC14-SN744_0143_BC01UKABXX.8_1_trim.fastq  /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/KIRC_data_and_RSEM_model/21a23647-5bce-4663-9c63-b1456c21503e/110715_UNC14-SN744_0143_BC01UKABXX.8_2_trim.fastq   /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Results/fusions/rsem/rsem-bowtie2/genv24_40event KIRC_21a23647_wRSPD 

# testing

~/Software/RSEM-1.2.28/rsem-simulate-reads /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Results/fusions/rsem/rsem-bowtie2/genv24_40event /work/DAT_086__TCGA_LAML/Staging/RNAseq/UNC_LCCC_data/21a23647-5bce-4663-9c63-b1456c21503e/KIRC_21a23647.stat/KIRC_21a23647.model /work/DAT_086__TCGA_LAML/Staging/RNAseq/UNC_LCCC_data/21a23647-5bce-4663-9c63-b1456c21503e/KIRC_21a23647.isoforms.results 029e-07 10000000 sim-test

# generate read data

module load rsem/1.2.28

rsem-simulate-reads /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Results/fusions/rsem/rsem-bowtie2/genv24_40event /work/DAT_086__TCGA_LAML/Staging/RNAseq/UNC_LCCC_data/21a23647-5bce-4663-9c63-b1456c21503e/KIRC_21a23647.stat/KIRC_21a23647.model /work/DAT_086__TCGA_LAML/Staging/RNAseq/UNC_LCCC_data/21a23647-5bce-4663-9c63-b1456c21503e/KIRC_21a23647.isoforms.results_modified 029e-07 50000000 sim1-50m


# align simulated data with star

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /external-data/Genome/indicies/Hsapiens_Gencode_GRCh38_STAR/ --genomeFastaFiles /external-data/Genome/genomes/Hsapiens_Gencode_GRCh38/primary/GRCh38.primary_assembly.genome.fa --sjdbGTFfile /external-data/Genome/gene_models/Hsapiens_gencode.v24.primary_assembly.annotation.gtf --sjdbOverhang 100

# assess alignments with STAR-fusion

DO NOT MODULE LOAD PERL

~/Software/STAR-Fusion_v0.6.0/STAR-Fusion --genome_lib_dir /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/GRCh38_gencode_v23_CTAT_lib/ -J /work/DAT_086__TCGA_LAML/Data/b3cf3b46-4585-4468-b2ed-15417b6d52e6Chimeric.out.junction --output_dir . --no_filter




# per-fusion event reference prep for rsem

module load bowtie rsem

rsem-prepare-reference --transcript-to-gene-map /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Results/fusions/rsem/gencode-fusion40-map-long.txt --bowtie2 --polyA  /external-data/Genome/transcriptomes/Hsapiens_Gencode_v24/gencode.v24.transcripts.fa,/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Results/fusions/FUSIM/2genes_exonboundaries_40events/fusions_2gene40exon.fasta genv24_40event



## creating diploid genome with phased SNPs

module load bcftools


# adjust chromosome names in phased SNP vcfs back to GRCh37
sed -e 's/^chr//g' filtered_germline_snv_HG002_hg19.0.vcf > filtered_germline_snv_HG002.0.vcf 

# compress and index vcf
bgzip filtered_germline_snv_HG002.0.vcf 
bcftools index filtered_germline_snv_HG002.0.vcf.gz 

# put variants in genome - hap1
bcftools consensus -f Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -o Homo_sapiens.GRCh37.75.dna.primary.hap1.fa ../../phased_SNP_SMC-het/filtered_germline_snv_HG002.0.vcf.gz

# Changing chromosome names in genome for diploid version
sed -e 's/>\([[:graph:]]*\)\s/>\1-hap1 /g' hap1/Homo_sapiens.GRCh37.75.dna.primary.hap1.fa > hap1/Homo_sapiens.GRCh37.75.dna.primary.hap1_wnames.fa

# Changing gene, transcript, and chromosome names in GTF
nohup ../../Scripts/rnaseq_fusion_simulation/make_haploid_gtf.py /external-data/Genome/gene_models/Hsapiens_Ensembl_v75.gtf hap1 > hap1/Hsapiens_Ensembl_v75_hap1.gtf &

# creating diploid reference
rsem-prepare-reference --gtf Hsapiens_Ensembl_v75_diploid.gtf --bowtie2 --polyA Homo_sapiens.GRCh37.75.primary.diploid.fa GRCH37v75_diploid &> prepref.log 


##########################
# generating sim1a dataset
##########################

java -jar ../../../../Scripts/fusim-0.2.2/fusim.jar --gene-model=/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/Hsapiens_Ensembl_v75.refflat.txt --fusions=10 --reference=/external-data/Genome/genomes/Hsapiens_Ensembl_GRCh37/primary_nomask/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --fasta-output=sim1a.fasta --text-output=sim1a_fusions.txt --cds-only --keep-exon-boundary &> sim1a_generation.log

# Get 10 transcripts that meet the length requirements and make a GTF
/Scripts/rnaseq_fusion_simulation/length_filter_fusim.py --number 10 unfiltered_sim1a.fasta > test.out

# Get separate fasta files for the 10 transcripts
/Scripts/rnaseq_fusion_simulation/parse_multi_fasta.py filtered_sim1a_fusions.fasta 

# Concatenate GTFs
cat /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/diploid_reference_genome/Hsapiens_Ensembl_v75_diploid.gtf unfiltered_sim1_filtered.gtf > Hsapiens_Ensembl_v75_pfusions.gtf 

# Make diploid genome + fusion reference
rsem-prepare-reference --gtf Hsapiens_Ensembl_v75_pfusions.gtf --star --num-threads 4 /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/diploid_reference_genome/Homo_sapiens.GRCh37.75.primary.diploid.fa,ENST00000260372-ENST00000373759.fasta,ENST00000300087-ENST00000525358.fasta,ENST00000355312-ENST00000453113.fasta,ENST00000393546-ENST00000538922.fasta,ENST00000442789-ENST00000432842.fasta,ENST00000457268-ENST00000567923.fasta,ENST00000548750-ENST00000446827.fasta,ENST00000556979-ENST00000411852.fasta,ENST00000563683-ENST00000530355.fasta,ENST00000565637-ENST00000452710.fasta GRCh37v75_diploid_sim1a &> prepref_sim1a.log 

# generate read data
module load rsem/1.2.28

rsem-simulate-reads GRCh37v75_diploid_sim1a /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/OICR_CPCG_prostate/star-ref/CPCG_0340.stat/CPCG_0340.model CPCG_0340.isoforms.results_modDiploid 0.066 30000000 sim1a_30m &>sim1a_30m.log

rsem-simulate-reads GRCh37v75_diploid_sim1a /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/OICR_CPCG_prostate/star-ref/CPCG_0340.stat/CPCG_0340.model CPCG_0340.isoforms.results_modDiploid_fusionsOnly 0.066 760 sim1a_30m_fus



## Generating separate references for sim1a

rsem-prepare-reference --gtf unfiltered_sim1_filtered.gtf --star --num-threads 4 ENST00000260372-ENST00000373759.fasta,ENST00000300087-ENST00000525358.fasta,ENST00000355312-ENST00000453113.fasta,ENST00000393546-ENST00000538922.fasta,ENST00000442789-ENST00000432842.fasta,ENST00000457268-ENST00000567923.fasta,ENST00000548750-ENST00000446827.fasta,ENST00000556979-ENST00000411852.fasta,ENST00000563683-ENST00000530355.fasta,ENST00000565637-ENST00000452710.fasta sim1a_10events &>prepref_sim1a_10events.log

rsem-prepare-reference --gtf ../Hsapiens_Ensembl_v75_diploid.gtf --star --num-threads 4 ../Homo_sapiens.GRCh37.75.primary.diploid.fa GRCh37v75_diploid &> prepref_diploid.log &


## Simulating reads

rsem-simulate-reads sim1a_10events /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/OICR_CPCG_prostate/star-ref/CPCG_0340.stat/CPCG_0340.model CPCG_0340.isoforms.results_modDiploid_fusionsOnly 0.066 790 sim1a_30m_fus

rsem-simulate-reads /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/diploid_reference_genome/STAR/GRCh37v75_diploid /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/OICR_CPCG_prostate/star-ref/CPCG_0340.stat/CPCG_0340.model CPCG_0340.isoforms.results_modDiploid_0 0.066 29000200 sim_diploid_30m &> sim_diploid_30m.log &

# merge read datasets
cat sim1a_30m_fus_2.fq sim_diploid_30m_2.fq | gzip > sim1a_30m_merged_2.fq.gz &