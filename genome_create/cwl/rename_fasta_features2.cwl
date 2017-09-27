#!/usr/bin/env cwl-runner
#
# Authors: Andrew Lamb

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [sed]

doc: "Changing chromosome names in genome for diploid version"

stdout: haploid2.fasta

inputs:

  input_string:
    type: string
    inputBinding:
      position: 1
      prefix: -e

  input_fasta:
    type: File
    inputBinding:
      position: 2

outputs:

  output_fasta:
    type: File
    outputBinding:
      glob: haploid2.fasta
