#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: Workflow

doc: "Fusion Simulation Workflow"

requirements:
  - class: MultipleInputFeatureRequirement

inputs:
  SIM_NAME: string
  GTF: File
  NUM_EVENTS: int
  TARGET_DEPTH: int
  GENOME: File
  EXPRESSION_PROFILE: File
  RSEM_MODEL: File
  DIP_GENOME: File
  SEED: ["null", int]
  MID_EXON_FUSIONS: ["null", boolean]

outputs:
  OUTPUT:
    type:
      type: array
      items: File
    outputSource: [fusion/fusionTruth, reads/isoformTruth, reads/fastq1, reads/fastq2, archive/archive]   
 
steps:

  tar:
    run: ../general_tools/tar_extract.cwl
    in:
      input: DIP_GENOME

    out: [output]

  gunzip:
    run: ../general_tools/gunzip.cwl
    in:
      input: GENOME

    out: [output]

  fusion:
    run: ../fusion_create/cwl/create_fusion.cwl
    in:
      gtf: GTF
      genome: gunzip/output
      numEvents: NUM_EVENTS
      simName: SIM_NAME
      seed: SEED
      mid_exon_fusions: MID_EXON_FUSIONS

    out: [fusGTF, fusRef, fusionTruth, fusLog, fusFA]

  isoform:
    run: ../model_isoforms/cwl/model_isoforms.cwl
    in:
      tpm: EXPRESSION_PROFILE
      gtf: fusion/fusGTF
      depth: TARGET_DEPTH
      seed: SEED

    out: [isoformTPM, fusionTPM, isoformLog]

  reads:
    run: ../fastq_create/cwl/create_fastq.cwl
    in:
      totalReads: TARGET_DEPTH
      simName: SIM_NAME
      RSEMmodel: RSEM_MODEL
      isoformTPM: isoform/isoformTPM
      fusionTPM: isoform/fusionTPM
      fusRef: fusion/fusRef
      dipGenome: tar/output
      isoformLog: isoform/isoformLog
      seed: SEED

    out: [isoformTruth, fastq1, fastq2, dip_gene_results, dip_iso_results, fus_gene_results, fus_iso_results, key1, key2] 

  archive:
    run: ../general_tools/tar_create.cwl
    in:
      fusion_log_file: fusion/fusLog
      fusion_FA_files: fusion/fusFA
      isoformTPM: isoform/isoformTPM
      fusionTPM: isoform/fusionTPM
      isoform_log_file: isoform/isoformLog
      dip_gene_results: reads/dip_gene_results
      dip_iso_results: reads/dip_iso_results
      fus_gene_results: reads/fus_gene_results
      fus_iso_results: reads/fus_iso_results
      key1: reads/key1
      key2: reads/key2



    out: [archive] 
