#!/usr/bin/env Rscript
# Originally written by Yin Hu in 2015 for Sage Bionetworks

setwd('/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/')

library(polyester)
library(Biostrings)
library(MASS)

## Set overall parameters
numreads = 5e7 # total number of reads to be generated
readLen = 100

nb_mu = 20 # mean of negative binomial
nb_theta = 2 # variance mu + mu^2/theta

## Read in transcript information
fastapath_fus = 'Data/fusion-simdata-rnaf/polyester_inputs/fusions.fasta'
numtx_fus = count_transcripts(fastapath_fus)
fastapath_ref = '/external-data/Genome/transcriptomes/Hsapiens_Gencode_v24/gencode.v24.pc_transcripts.fa'
numtx_ref = count_transcripts(fastapath_ref)
numtx = numtx_fus + numtx_ref
translen_fus = read.delim('Data/fusion-simdata-rnaf/polyester_inputs/translen_fusion.txt', header=F); colnames(translen_fus) = c("ID", "Length")
translen_ref = fasta.seqlengths(fastapath_ref,use.names = FALSE)


## Randomly assign transcript abundance
abundance = rnegbin(numtx, mu = nb_mu, theta = nb_theta) + 1
abundance = abundance * c(translen_fus$Length, translen_ref) # total coverage
abundance = abundance / sum(abundance)
reads = ceiling(numreads * abundance)
reads_fus = reads[1:numtx_fus]
reads_ref = reads[(1+numtx_fus):numtx]

#### Directly use the fusion fasta file
file.create("Results/fusions/simulated_reads_ref_only/finalsample_01_1.fasta")
file.create("Results/fusions/simulated_reads_ref_only/finalsample_01_2.fasta")
#file.create("Results/fusions/simulated_reads_ref_only/finalsample_02_1.fasta")
#file.create("Results/fusions/simulated_reads_ref_only/finalsample_02_2.fasta")
gc()
numshard = 20
# This package is not scalable. It will crash if simulating too many reads. Separate into multiple runs (10 here), and merge the results. Flaw here: sum of all shards may exceed the total number of reads for each transcript, due to the ceiling calc.
for (i in 1:numshard) { 
#  simulate_experiment(fastapath_ref, reads_per_transcript=ceiling(reads_ref/numshard), num_reps=1,fold_changes=rep(1, numtx_ref), outdir='Results/fusions/simulated_reads_ref_only/', transcriptid=rep("ref", numtx_ref), seed=12) # simple
  
  simulate_experiment(fastapath_ref, reads_per_transcript=ceiling(reads_ref/numshard), num_reps=1,fold_changes=rep(1, numtx_ref), outdir='Results/fusions/simulated_reads_ref_only/', transcriptid=rep("ref", numtx_ref), seed=12, distr='empirical',error_model='illumina5',bias='cdnaf',readlen=readLen) # noisier
  
  file.append("Results/fusions/simulated_reads_ref_only/finalsample_01_1.fasta", "Results/fusions/simulated_reads_ref_only/sample_01_1.fasta")
  file.append("Results/fusions/simulated_reads_ref_only/finalsample_01_2.fasta", "Results/fusions/simulated_reads_ref_only/sample_01_2.fasta")
}
gc()
simulate_experiment(fastapath_fus, reads_per_transcript=reads_fus, num_reps=1, fold_changes=rep(1, numtx_fus), outdir='Results/fusions/simulated_reads_fusion_only/', transcriptid=rep("syn", numtx_fus), seed=12)

