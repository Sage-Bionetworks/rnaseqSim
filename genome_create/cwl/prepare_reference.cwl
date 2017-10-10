#!/usr/bin/env cwl-runner
#
# Authors: Andrew Lamb

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [rsem-prepare-reference, --star]

doc: "Create diploid reference files"

hints:
  DockerRequirement:
    dockerPull: andrewlambsage/rnaseqsim:gw

requirements:
  - class: InlineJavascriptRequirement

inputs:

  input_gtf:
    type: File
    inputBinding:
      position: 1
      prefix: --gtf

  num_cores:
    type: File
    inputBinding:
      position: 2
      prefix: --num-threads

  input_fasta:
    type: File
    inputBinding:
      position: 3
  
  ref_name:
    type: string
    inputBinding:
      position: 4

outputs:

  output_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: $(inputs.ref_name + "*")
