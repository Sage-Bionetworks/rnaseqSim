#!/bin/bash


SIM_NAME=$1
NUM_EVENTS=$2
DEPTH=$3
EXPRESION_PROFILE=$4
RSEM_MODEL=$5
TIMEOUT=$6

CWL_PATH=rnaseq_fusion_simulation/workflow/fusion_simulation_workflow.cwl
INPUT_JOB=$SIM_NAME.json


./rnaseq_fusion_simulation/generate_job.py --syn-table SMC-RNA-Eval/syn.table $CONTEST_SIGN $CWL_PATH $TUMOR_ID > $INPUT_JOB
cwltool $CWL_PATH $INPUT_JOB

if [ "$6" != "" ]; then
  sudo poweroff
fi
