#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

doc: "Diploid Genome Creation Workflow"

requirements:
  - class: MultipleInputFeatureRequirement

inputs:
  INPUT_HAPLOID_FILE: File
  HAPOLID_COPY_NAME: string
  INPUT_VCF_FILE1: File
  INPUT_VCF_FILE2: File
  VARIANT_HAPLOID_FILE1: string
  VARIANT_HAPLOID_FILE2: string
  SED_STRING1: string
  SED_STRING2: string
  INPUT_GTF_FILE: File
  GTF_SUFFIX_STRING1: string
  GTF_SUFFIX_STRING2: string
  REF_NAME: string

outputs:
  OUTPUT:
    type: File
    outputSource: tar_create/archive
 
steps:

  copy_genome:
    run: copy_genome.cwl
    in: 
      input_fasta: INPUT_HAPLOID_FILE
      output_fasta_name: HAPOLID_COPY_NAME
    out: [output_fasta]

  index_vcf1:
    run: index_vcf.cwl
    in: 
      input_vcf: INPUT_VCF_FILE1
    out: [output_index]

  index_vcf2:
    run: index_vcf.cwl
    in: 
      input_vcf: INPUT_VCF_FILE2
    out: [output_index]

  add_variants1:
    run: add_variants.cwl
    in: 
      input_fasta: INPUT_HAPLOID_FILE
      output_fasta_name: VARIANT_HAPLOID_FILE1
      input_vcf: index_vcf1/output_index
    out: [output_fasta]

  add_variants2:
    run: add_variants.cwl
    in: 
      input_fasta: copy_genome/output_fasta
      output_fasta_name: VARIANT_HAPLOID_FILE2
      input_vcf: index_vcf2/output_index
    out: [output_fasta]

  rename_fasta_features1:
    run: rename_fasta_features1.cwl
    in: 
      input_string: SED_STRING1
      input_fasta: add_variants1/output_fasta
    out: [output_fasta]

  rename_fasta_features2:
    run: rename_fasta_features2.cwl
    in: 
      input_string: SED_STRING2
      input_fasta: add_variants2/output_fasta
    out: [output_fasta]

  combine_fastas:
    run: combine_fastas.cwl
    in: 
      input_fasta1: rename_fasta_features1/output_fasta
      input_fasta2: rename_fasta_features2/output_fasta
    out: [output_fasta]

  rename_gtf_features1:
    run: rename_gtf_features1.cwl
    in: 
      input_string: GTF_SUFFIX_STRING1
      input_gtf: INPUT_GTF_FILE
    out: [output_gtf]

  rename_gtf_features2:
    run: rename_gtf_features2.cwl
    in: 
      input_string: GTF_SUFFIX_STRING2
      input_gtf: INPUT_GTF_FILE
    out: [output_gtf]

  combine_gtfs:
    run: combine_gtfs.cwl
    in: 
      input_gtf1: rename_gtf_features1/output_gtf
      input_gtf2: rename_gtf_features2/output_gtf
    out: [output_gtf]
  
  prepare_reference:
    run: prepare_reference.cwl
    in: 
      input_gtf: combine_gtfs/output_gtf
      input_fasta: combine_fastas/output_fasta
      ref_name: REF_NAME
    out: [output_files]

  tar_create:
    run: tar_create.cwl
    in: 
      files: prepare_reference/output_files
    out: [archive]
