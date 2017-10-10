#!/usr/bin/env cwl-runner
#
# Authors: Andrew Lamb

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [tar, cvzf, archive.tgz]

doc: "command line: tar"

inputs:

  fusion_log_file:
    type: File
    inputBinding:
      position: 1

  fusion_FA_files:
    type:
      type: array
      items: File
    inputBinding:
      position: 2

  isoformTPM:
    type: File
    inputBinding:
      position: 3

  fusionTPM:
    type: File
    inputBinding:
      position: 4

  isoform_log_file:
    type: File
    inputBinding:
      position: 5

  dip_gene_results:
    type: File
    inputBinding:
      position: 6

  dip_iso_results:
    type: File
    inputBinding:
      position: 7

  fus_gene_results:
    type: File
    inputBinding:
      position: 8

  dip_iso_results:
    type: File
    inputBinding:
      position: 9

  fus_iso_results:
    type: File
    inputBinding:
      position: 10

  key1:
    type: File
    inputBinding:
      position: 11

  key2:
    type: File
    inputBinding:
      position: 12

outputs:

  archive:
    type: File
    outputBinding:
      glob: archive.tgz
