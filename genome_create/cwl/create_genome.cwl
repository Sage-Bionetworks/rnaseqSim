#!/usr/bin/env cwl-runner
#
# Authors: Allison Creason, Kristen Dang, Kyle Ellrott, Ryan Spangler

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [make_haploid_gtf.py]
stdout: haploid.gtf 

doc: "Create genome"

hints:
  DockerRequirement:
    dockerPull: alliecreason/create_genome

requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 80000

inputs:

  gtf:
    type: File
    inputBinding:
      position: 1

  suffix:
    type: string
    inputBinding:
      position: 2

outputs:

  output:
    type: stdout

