#!/usr/bin/env cwl-runner
#
# Authors: Andrew Lamb

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [tar, cvzf, archive.tgz]

doc: "command line: tar"

inputs:

  fusion_gtf_file:
    type: File
    inputBinding:
      position: 1

  fusion_ref_file:
    type: File
    inputBinding:
      position: 2

outputs:

  archive:
    type: File
    outputBinding:
      glob: archive.tgz
