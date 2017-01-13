# rnaseqSim
For SMC-RNA challenge, code for analyzing datasets of validated fusions and generating synthetic fusion data.

## Requirements:

STAR 2.4.2a (https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz)

RSEM v1.2.31 (https://github.com/deweylab/RSEM/archive/v1.2.31.tar.gz)

## Steps:

Step 1: Create fusion events, truth file, and RSEM-format fusion reference: `fusion_create/module.py`

Step 2: Adjust estimated isoform values to include expression for fusion genes according to a model: `model_isoforms/modify_model_tpm_for_diploid.R`

Step 3: Generate reads from diploid and fusion references: `fastq_create/generate_reads_noToil.py`


## To run cwl workflow:

`cwltool workflow/fusion_simulation_workflow.cwl [INPUT.JSON]`

