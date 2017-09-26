#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [/home/aelamb/Projects/create_genome/misc/prepare_files.sh]

doc: "Create diploid reference files"

inputs:
 
  input_string:
    type: string
    inputBinding:
      position: 1

outputs:

  output_array:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.txt"

