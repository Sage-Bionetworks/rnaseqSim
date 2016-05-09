#! /bin/bash
#$ -N rsem
#$ -pe threads 8
#$ -q regular
#$ -S /bin/bash
#$ -l h_vmem=32G
#$ -j y
#$ -V
#$ -cwd


module load rsem star/2.4.2a

sample=`basename $1`

rsem-calculate-expression --num-threads 8 --paired-end --estimate-rspd --strand-specific --star $1 $2 /external-data/Genome/indicies/Hsapiens_Ensembl_GRCh37_RSEM/STAR_ENSGv75/GRCh37v75_STAR $sample &> $sample.align.log

