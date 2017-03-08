#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: Workflow

doc: "Fusion Simulation Workflow"

inputs:
  SIM_NAME: string
  GTF: File
  NUM_EVENTS: int
  TARGET_DEPTH: int
  GENOME: File
  EXPRESSION_PROFILE: File
  RSEM_MODEL: File
  DIP_GENOME: File

outputs:
  OUTPUT:
    type: File
    outputSource: [reads/fusionTruth, reads/isoformTruth, reads/fastq1, reads/fastq2]

steps:

  genome:
    run: ../genome_create/cwl/create_genome.cwl
    in:
      gtf: GTF
      suffix: SIM_NAME

    out: [output]

  fusion:
    run: ../fusion_create/cwl/create_fusion.cwl
    in:
      gtf: GTF
      genome: GENOME
      numEvents: NUM_EVENTS
      simName: SIM_NAME

    out: [fusGTF, fusRef]

  isoform:
    run: ../model_isoforms/cwl/model_isoforms.cwl
    in:
      tpm: EXPRESSION_PROFILE
      gtf: fusion/fusGTF
      depth: TARGET_DEPTH

    out: [isoformTPM, fusionTPM, log]

  reads:
    run: ../fastq_create/cwl/create_fastq.cwl
    in:
      totalReads: TARGET_DEPTH
      simName: SIM_NAME
      RSEMmodel: RSEM_MODEL
      isoformTPM: isoform/isoformTPM
      fusionTPM: isoform/fusionTPM
      fusRef: fusion/fusRef
      dipGenome: DIP_GENOME
      isoformLog: isoform/log

    out: [fusionTruth, isoformTruth, fastq1, fastq2]
