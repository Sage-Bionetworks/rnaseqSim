#!/usr/bin/env cwl-runner
#
# Authors: Allison Creason, Kristen Dang, Kyle Ellrott, Ryan Spangler

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [generate_reads_noToil.py]

doc: "Generate FastQ reads files"

hints:
  DockerRequirement:
    dockerPull: alliecreason/create_fastq

requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 
    ramMin: 

inputs:

  totalReads:
    type: int?
    inputBinding:
      position: 1
      prefix: --totalReads

  numSimReads:
    type: int?
    inputBinding:
      position: 1
      prefix: --numSimReads
  
  simName:
    type: string
    inputBinding:
      position: 1
      prefix: --simName

  RSEMmodel:
    type: File
    inputBinding:
      position: 1
      prefix: --RSEMmodel

  isoformTPM:
    type: File
    inputBinding:
      position: 1
      prefix: --isoformTPM

  fusionTPM:
    type: File
    inputBinding:
      position: 1
      prefix: --fusionTPM

outputs:

  output:
    type: File

arguments:
  - valueFrom: $(inputs.simName + "_" + inputs.numEvents +  "_ev")
    prefix: "--fusRef"
