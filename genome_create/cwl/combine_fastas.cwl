#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [cat]

doc: "Combine fasta files"

stdout: diploid.fasta

inputs:

  input_fasta1:
    type: File
    inputBinding:
      position: 1

  input_fasta2:
    type: File
    inputBinding:
      position: 2

outputs:

  output_fasta:
    type: File
    outputBinding:
      glob: diploid.fasta
