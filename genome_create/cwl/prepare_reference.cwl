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

inputs:

  input_gtf:
    type: File
    inputBinding:
      position: 1
      prefix: --gtf

  input_fasta:
    type: File
    inputBinding:
      position: 2
  
  ref_name:
    type: string
    inputBinding:
      position: 3

outputs:

  output_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: $(inputs.ref_name + "*")
