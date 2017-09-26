#!/usr/bin/env cwl-runner
#
# Authors: Andrew Lamb

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [tar, cvzf, genome.tgz]

doc: "command line: tar"

inputs:
 
  files:
    type:
      type: array
      items: File
      inputBinding:
        separate: true
    inputBinding:
      position: 1
 
outputs:

  archive:
    type: File
    outputBinding:
      glob: genome.tgz
