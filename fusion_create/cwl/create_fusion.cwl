#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [module.py]

doc: "Create fusion events"

hints:
  DockerRequirement:
    dockerPull: alliecreason/rnaseqsim

requirements:
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 
    ramMin: 

stdout: fus.log

inputs:

  gtf:
    type: File
    inputBinding:
      position: 1
      prefix: --gtf

  genome:
    type: File
    inputBinding:
      position: 1
      prefix: --genome

  numEvents:
    type: int
    inputBinding:
      position: 1
      prefix: --numEvents

  simName:
    type: string
    inputBinding:
      position: 1
      prefix: --simName
  
  seed:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --seed

  mid_exon_fusions:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: --mid_exon_fusions

outputs:

  fusionTruth:
    type: File
    outputBinding:
      glob: $(inputs.simName + "_filtered.bedpe")

  fusGTF:
    type: File
    outputBinding:
      glob: $(inputs.simName + ".gtf")

  fusLog:
    type: stdout

  fusFA:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.fasta"

  fusRef:
    type: File
    outputBinding:
      glob: $(inputs.simName + "_" + inputs.numEvents +  "_ev.seq")
    secondaryFiles:
     - ^.chrlist
     - ^.grp
     - ^.idx.fa
     - ^.n2g.idx.fa
     - ^.ti
     - ^.transcripts.fa
