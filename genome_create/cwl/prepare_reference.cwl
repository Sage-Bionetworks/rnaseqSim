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
    type: int
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

  output_chrlist:
    type: File
    outputBinding:
      glob: $(inputs.ref_name + ".chrlist")

  output_grp:
    type: File
    outputBinding:
      glob: $(inputs.ref_name + ".grp")

  output_Log.out:
    type: File
    outputBinding:
      glob: $(inputs.ref_name + "Log.out")

  output_n2g.idx.fa:
    type: File
    outputBinding:
      glob: $(inputs.ref_name + ".n2g.idx.fa")

  output_idx.fa:
    type: File
    outputBinding:
      glob: $(inputs.ref_name + ".idx.fa")

  output_seq:
    type: File
    outputBinding:
      glob: $(inputs.ref_name + ".seq")

  output_ti:
    type: File
    outputBinding:
      glob: $(inputs.ref_name + ".ti")

  output_transcripts.fa:
    type: File
    outputBinding:
      glob: $(inputs.ref_name + ".transcripts.fa")

