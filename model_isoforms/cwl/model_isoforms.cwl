#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [modify_model_tpm_for_diploid.R]

doc: "Model Isoform Expression"

hints:
  DockerRequirement:
    dockerPull: alliecreason/rnaseqsim

requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 
    ramMin: 

stdout: isoform.log

inputs:

  tpm:
    type: File
    inputBinding:
      position: 1
      prefix: --TPM

  gtf:
    type: File
    inputBinding:
      position: 1
      prefix: --gtf
  
  depth:
    type: int
    inputBinding:
      position: 1
      prefix: --targetDepth

outputs:

  isoformTPM:
    type: File
    outputBinding:
      glob: $("*results_modDiploid_" + inputs.depth)

  fusionTPM:
    type: File
    outputBinding:
      glob: $("*results_modDiploidFusionOnly_" + inputs.depth)

  isoformLog:
    type: stdout
