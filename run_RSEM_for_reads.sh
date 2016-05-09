#! /bin/bash
#$ -N rsem
#$ -q regular
#$ -S /bin/bash
#$ -l h_vmem=15G
#$ -j y
#$ -V
#$ -cwd


#$1 = RSEM model file
#$2 = non-fusion isoform values
#$3 = number of non-fusion reads
#$4 = simulation name
#$5 = fusion reference
#$6 = fusion isoform values
#$7 = number of fusion reads

module load rsem/1.2.8 star/2.4.2a



rsem-simulate-reads /work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data/diploid_reference_genome/STAR/GRCh37v75_diploid $1 $2 0.066 $3 $4.diploid 

rsem-simulate-reads $5 $1 $6 0.066 $7 $4.fusions 

# merge read datasets
cat $4.diploid_1.fq $4.fusions_1.fq | gzip > $4_merged_1.fq.gz 
cat $4.diploid_2.fq $4.fusions_2.fq | gzip > $4_merged_2.fq.gz 
