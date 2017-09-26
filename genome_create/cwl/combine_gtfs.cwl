#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [cat]

doc: "Combine gtf files"

stdout: diploid.gtf

inputs:

  input_gtf1:
    type: File
    inputBinding:
      position: 1

  input_gtf2:
    type: File
    inputBinding:
      position: 2

outputs:

  output_gtf:
    type: File
    outputBinding:
      glob: diploid.gtf
