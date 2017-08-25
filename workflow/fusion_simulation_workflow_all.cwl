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

outputs:
  OUTPUT:
    type:
      type: array
      items: File
    outputSource: [fusion/fusionTruth, reads/isoformTruth, reads/fastq1, reads/fastq2, archive/archive]   
 
steps:

  tar:
    run: ../general_tools/tar.cwl
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

    out: [fusGTF, fusRef, fusionTruth]

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

    out: [isoformTruth, fastq1, fastq2] 

  archive:
    run: ../general_tools/tar2.cwl
    in:
      fusion_gtf_file: fusion/fusGTF
      fusion_ref_file: fusion/fusRef
      
    out: [archive] 


