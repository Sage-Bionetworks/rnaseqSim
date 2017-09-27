#!/bin/bash
wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
wget http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
wget http://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20170905.vcf.gz
cp clinvar_20170905.vcf.gz clinvar_20170905.vcf2.gz
