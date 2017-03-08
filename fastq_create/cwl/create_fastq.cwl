#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [generate_reads.py]

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
      valueFrom: $(inputs.fusRef + "/" + inputs.fusRef.nameroot)

  dipGenome:
    type: File
    inputBinding:
      position: 1
      prefix: --dipGenome

  isoformLog:
    type: File
    inputBinding:
      position: 1
      prefix: --isoformLog

outputs:

  fastq1:
    type: File
    outputBinding:
      glob: $(simName + "_mergeSort_1.fq.gz")

  fastq2:
    type: File
    outputBinding:
      glob: $(simName + "_mergeSort_2.fq.gz")

  isoformTruth:
    type: File
    outputBinding:
      glob: $(simName + "_isoforms_truth.txt")

  fusionTruth:
    type: File
    outputBinding:
      glob: $(simName + "_filtered.bedpe")
