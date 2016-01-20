#!/usr/bin/env Rscript
# Originally written by Yin Hu in 2015 for Sage Bionetworks


library(polyester)
library(Biostrings)
library(MASS)

## Set overall parameters
numreads = 50000000 # total number of reads to be generated

nb_mu = 20 # mean of negative binomial
nb_theta = 2 # variance mu + mu^2/theta

## Read in transcript information
fastapath_fus = "Fusim//fusions.fasta"
numtx_fus = count_transcripts(fastapath_fus)
fastapath_ref = "ref.fa"
numtx_ref = count_transcripts(fastapath_ref)
numtx = numtx_fus + numtx_ref
translen_fus = read.delim("translen_fusion.txt", header=F); colnames(translen_fus) = c("ID", "Length")
translen_ref = read.delim("translen_ref.txt", header=F); colnames(translen_ref) = c("ID", "Length")

## Randomly assign transcript abundance
abundance = rnegbin(numtx, mu = nb_mu, theta = nb_theta) + 1
abundance = abundance * c(translen_fus$Length, translen_ref$Length) # total coverage
abundance = abundance / sum(abundance)
reads = ceiling(numreads * abundance)
reads_fus = reads[1:numtx_fus]
reads_ref = reads[(1+numtx_fus):numtx]

#### Directly use the fusion fasta file
file.create("simulated_reads_ref_only/finalsample_01_1.fasta")
file.create("simulated_reads_ref_only/finalsample_01_2.fasta")
file.create("simulated_reads_ref_only/finalsample_02_1.fasta")
file.create("simulated_reads_ref_only/finalsample_02_2.fasta")
numshard = 10
for (i in 1:numshard) { # This package is not scalable. It will crash if simulating too many reads. Separate into multiple runs (10 here), and merge the results
  simulate_experiment(fastapath_ref, reads_per_transcript=ceiling(reads_ref/numshard), num_reps=1,
                      fold_changes=rep(1, numtx_ref), outdir='simulated_reads_ref_only/', transcriptid=rep("ref", numtx_ref), seed=12) # Flaw here: sum of all shards may exceed the total number of reads for each transcript, due to the ceiling calc
  file.append("simulated_reads_ref_only/finalsample_01_1.fasta", "simulated_reads_ref_only/sample_01_1.fasta")
  file.append("simulated_reads_ref_only/finalsample_01_2.fasta", "simulated_reads_ref_only/sample_01_2.fasta")
  file.append("simulated_reads_ref_only/finalsample_02_1.fasta", "simulated_reads_ref_only/sample_02_1.fasta")
  file.append("simulated_reads_ref_only/finalsample_02_2.fasta", "simulated_reads_ref_only/sample_02_2.fasta")
}

simulate_experiment(fastapath_fus, reads_per_transcript=reads_fus, num_reps=1,
                    fold_changes=rep(1, numtx_fus), outdir='simulated_reads_fusion_only/', transcriptid=rep("syn", numtx_fus), seed=12)

