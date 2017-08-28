# rnaseqSim
For SMC-RNA challenge, code for analyzing datasets of validated fusions and generating synthetic fusion data.

## CWL Workflow:

Running the workflow requires CWL v1.0+ and Docker.

To run the workflow:

`cwltool workflow/fusion_simulation_workflow.cwl [INPUT.JSON]`
`cwltool workflow/fusion_simulation_workflow_all.cwl [INPUT.JSON]`

The input JSON needs the fields:

  SIM_NAME: string
  GTF: File
  NUM_EVENTS: int
  TARGET_DEPTH: int
  GENOME: File
  EXPRESSION_PROFILE: File
  RSEM_MODEL: File
  DIP_GENOME: File

And can use the optional fields
  SEED: ["null", int]


## Description of inputs

SIM_NAME:
GTF:
NUM_EVENTS:
TARGET_DEPTH:
GENOME:
EXPRESSION_PROFILE:
RSEM_MODEL:
DIP_GENOME:

SEED: This is optional. If given all scripts with a random element in the 
workflow will have a seed set at the given integer.


## Description of outputs

[SIM_NAME]_filtered.bedpe:

[SIM_NAME]_isoforms_truth.txt:

[SIM_NAME]_mergeSort_1.fq.gz:

[SIM_NAME]_mergeSort_2.fq.gz:

archive.tgz: This will store other intermeduate files if 
`fusion_simulation_workflow_all.cwl`was used.


# older descriptions:

## Requirements:

[STAR 2.4.2a] (https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz)

[RSEM v1.2.31] (https://github.com/deweylab/RSEM/archive/v1.2.31.tar.gz)

## Required Inputs

Diploid Genome - Homo_sapiens.GRCh37.75.primary.diploid.fa.gz (syn8348583)

Diploid GTF - Hsapiens_Ensembl_v75_diploid.gtf.gz (syn8348617)

Reference GTF - Hsapiens_Ensembl_v75_refonly.gtf (syn8348668)

Model file - CPCG_0258.R1.fastq.model (syn8348382)

Expression profile -


## Basic Steps:

Step 1 - Index Diploid Genome:

`rsem-prepare-reference --gtf [diploid.ref.gtf] --star [diploid.ref.fa] [Index name]`

Step 2 - Create fusion events, truth file, and RSEM-format fusion reference:

`fusion_create/module.py --gtf Hsapiens_Ensembl_v75_refonly.gtf --numEvents [XX] --simName [simName]`

Step 3 - Adjust estimated isoform values to include expression for fusion genes according to a model:

`model_isoforms/modify_model_tpm_for_diploid.R --TPM [input expression profile] --gtf [simName.gtf] --targetDepth [XX] --codeDir [/path/to/code] &> [output.log]`

Step 4 - Generate reads from diploid and fusion references:

`fastq_create/generate_reads.py --totalReads [targetDepth * 1000000] --numSimReads [output.log] --simName [simName] --RSEMmodel [model file] --isoformTPM [model_isoforms output] --fusionTPM [model_isoforms output] --fusRef [fusion_create output]`

