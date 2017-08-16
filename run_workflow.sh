#!/bin/bash

SIM_NAME=$1
NUM_EVENTS=$2
DEPTH=$3
EXPRESION_PROFILE=$4
RSEM_MODEL=$5
TIMEOUT=$6

HOME="/home/ubuntu"
INPUT_PATH="$HOME/inputs"
OUTPUT_PATH="$HOME/outputs"
RNASEQSIM_PATH="$HOME/rnaseqSim"
OUTPUT_BUCKET="gs://smc-rna-eval/training/"

if [ ! -d "$INPUT_PATH" ]; then
       mkdir $INPUT_PATH
fi

INPUT_JOB="$INPUT_PATH/input.json"
$RNASEQSIM_PATH/generate_job.py --simName $SIM_NAME --numEvents $NUM_EVENTS --targetDepth $DEPTH --expressionProfile $EXPRESION_PROFILE --RSEMmodel $RSEM_MODEL 

if [ ! -d "$OUTPUT_PATH" ]; then
       mkdir $OUTPUT_PATH
fi

CWL_PATH="$RNASEQSIM_PATH/workflow/fusion_simulation_workflow.cwl"
cd $OUTPUT_PATH
cwltool $CWL_PATH $INPUT_JOB
gsutil cp $OUTPUT_PATH/* $OUTPUT_BUCKET


#if [ "$6" != "" ]; then
#  sudo poweroff
#fi
