#!/usr/bin/env cwl-runner
#
# Authors: Andrew Lamb

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [tar, cvzfh, genome.tgz]

doc: "command line: tar"

inputs:

  input_chrlist:
    type: File
    inputBinding:
      position: 1

  input_grp:
    type: File
    inputBinding:
      position: 2

  input_Log.out:
    type: File
    inputBinding:
      position: 3

  input_n2g.idx.fa:
    type: File
    inputBinding:
      position: 4

  input_idx.fa:
    type: File
    inputBinding:
      position: 5

  input_seq:
    type: File
    inputBinding:
      position: 6

  input_ti:
    type: File
    inputBinding:
      position: 7

  input_transcripts.fa:
    type: File
    inputBinding:
      position: 8
 
outputs:

  archive:
    type: File
    outputBinding:
      glob: genome.tgz
