#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [generate_reads.py]

doc: "Generate FastQ reads files"

hints:
  DockerRequirement:
    dockerPull: andrewelambsage/rnaseqsim

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
      valueFrom: $(inputs.fusRef.dirname + "/" + inputs.fusRef.nameroot)

  dipGenome:
    type: Directory
    inputBinding:
      position: 1
      prefix: --dipGenome
      valueFrom: $(inputs.dipGenome.path + "/GRCh37v75_STAR/GRCh37v75_STAR")
      
  isoformLog:
    type: File
    inputBinding:
      position: 1
      prefix: --isoformLog
      
  seed:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --seed

outputs:

  fastq1:
    type: File
    outputBinding:
      glob: $(inputs.simName + "_mergeSort_1.fq.gz") 

  fastq2:
    type: File
    outputBinding:
      glob: $(inputs.simName + "_mergeSort_2.fq.gz") 

  isoformTruth:
    type: File
    outputBinding:
      glob: $(inputs.simName + "_isoforms_truth.txt") 

