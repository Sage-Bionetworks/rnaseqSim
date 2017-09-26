#!/usr/bin/env cwl-runner
#
# Authors: Andrew Lamb

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [bcftools, index]

doc: "Index VCF file"

hints:
  DockerRequirement:
    dockerPull: rnaseqsim

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.input_vcf)

inputs:

  input_vcf:
    type: File
    inputBinding:
      position: 1

outputs:

  output_index:
    type: File
    outputBinding:
      glob: $(inputs.input_vcf.basename)
    secondaryFiles:
      - .csi

arguments:
  - valueFrom: $(inputs.input_vcf.basename + ".csi")
    position: 2
  


