#!/usr/bin/env cwl-runner
#
# Authors: Andrew Lamb

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [make_haploid_gtf.py]

doc: "Changing gene, transcript, and chromosome names in GTF"

hints:
  DockerRequirement:
    dockerPull: andrewlambsage/rnaseqsim:gw

hints:
  DockerRequirement:
    dockerPull: rnaseqsim

stdout: haploid1.gtf

inputs:

  input_string:
    type: string
    inputBinding:
      position: 1
      prefix: --suffix

  input_gtf:
    type: File
    inputBinding:
      position: 2
      prefix: --GTF

outputs:

  output_gtf:
    type: File
    outputBinding:
      glob: haploid1.gtf
