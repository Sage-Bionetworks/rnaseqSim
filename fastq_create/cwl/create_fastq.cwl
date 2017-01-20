#!/usr/bin/env cwl-runner

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

  totalReads:
    type: int
    inputBinding:
      position: 1
      prefix: --totalReads
      valueFrom: $(inputs.totalReads * 1000000)

  fusRef:
    type: File
    inputBinding:
      position: 1
      prefix: --fusRef
      valueFrom: $(inputs.fusRef.nameroot)

  dipGenome:
    type: File
    inputBinding:
      position: 1
      prefix: --dipGenome

outputs:

  output:
    type: File

