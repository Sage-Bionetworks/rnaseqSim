#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [cp]

doc: "Make copy of genome"

inputs:

  input_fasta:
    type: File
    inputBinding:
      position: 1

  output_fasta_name:
    type: string
    inputBinding:
      position: 2

outputs:

  output_fasta:
    type: File
    outputBinding:
      glob: $(inputs.output_fasta_name)


