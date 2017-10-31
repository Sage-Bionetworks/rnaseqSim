#!/usr/bin/env cwl-runner
#
# Authors: Andrew Lamb

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
  NUM_CORES: int

outputs:
  OUTPUT:
    type:
      type: array
      items: File
    outputSource: [tar_create/archive, combine_gtfs/output_gtf, combine_fastas/output_fasta]  
 
steps:

  copy_genome:
    run: ../genome_create/cwl/copy_genome.cwl
    in: 
      input_fasta: INPUT_HAPLOID_FILE
      output_fasta_name: HAPOLID_COPY_NAME
    out: [output_fasta]

  index_vcf1:
    run: ../genome_create/cwl/index_vcf.cwl
    in: 
      input_vcf: INPUT_VCF_FILE1
    out: [output_index]

  index_vcf2:
    run: ../genome_create/cwl/index_vcf.cwl
    in: 
      input_vcf: INPUT_VCF_FILE2
    out: [output_index]

  add_variants1:
    run: ../genome_create/cwl/add_variants.cwl
    in: 
      input_fasta: INPUT_HAPLOID_FILE
      output_fasta_name: VARIANT_HAPLOID_FILE1
      input_vcf: index_vcf1/output_index
    out: [output_fasta]

  add_variants2:
    run: ../genome_create/cwl/add_variants.cwl
    in: 
      input_fasta: copy_genome/output_fasta
      output_fasta_name: VARIANT_HAPLOID_FILE2
      input_vcf: index_vcf2/output_index
    out: [output_fasta]

  rename_fasta_features1:
    run: ../genome_create/cwl/rename_fasta_features1.cwl
    in: 
      input_string: SED_STRING1
      input_fasta: add_variants1/output_fasta
    out: [output_fasta]

  rename_fasta_features2:
    run: ../genome_create/cwl/rename_fasta_features2.cwl
    in: 
      input_string: SED_STRING2
      input_fasta: add_variants2/output_fasta
    out: [output_fasta]

  combine_fastas:
    run: ../genome_create/cwl/combine_fastas.cwl
    in: 
      input_fasta1: rename_fasta_features1/output_fasta
      input_fasta2: rename_fasta_features2/output_fasta
    out: [output_fasta]

  rename_gtf_features1:
    run: ../genome_create/cwl/rename_gtf_features1.cwl
    in: 
      input_string: GTF_SUFFIX_STRING1
      input_gtf: INPUT_GTF_FILE
    out: [output_gtf]

  rename_gtf_features2:
    run: ../genome_create/cwl/rename_gtf_features2.cwl
    in: 
      input_string: GTF_SUFFIX_STRING2
      input_gtf: INPUT_GTF_FILE
    out: [output_gtf]

  combine_gtfs:
    run: ../genome_create/cwl/combine_gtfs.cwl
    in: 
      input_gtf1: rename_gtf_features1/output_gtf
      input_gtf2: rename_gtf_features2/output_gtf
    out: [output_gtf]
  
  prepare_reference:
    run: ../genome_create/cwl/prepare_reference.cwl
    in: 
      input_gtf: combine_gtfs/output_gtf
      num_cores: NUM_CORES
      input_fasta: combine_fastas/output_fasta
      ref_name: REF_NAME
    out: [output_chrlist, output_grp, output_Log.out, output_n2g.idx.fa, output_idx.fa, output_seq, output_ti, output_transcripts.fa]

  tar_create:
    run: ../genome_create/cwl/tar_create.cwl
    in: 
      input_chrlist: prepare_reference/output_chrlist
      input_grp: prepare_reference/output_grp
      input_Log.out: prepare_reference/output_Log.out
      input_n2g.idx.fa: prepare_reference/output_n2g.idx.fa
      input_idx.fa: prepare_reference/output_idx.fa
      input_seq: prepare_reference/output_seq
      input_ti: prepare_reference/output_ti
      input_transcripts.fa: prepare_reference/output_transcripts.fa
    out: [archive]
