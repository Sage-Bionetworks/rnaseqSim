#!/usr/bin/env python

import os
import yaml
import argparse
import json
import synapseclient
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--simName")
    parser.add_argument("--numEvents")
    parser.add_argument("--targetDepth")
    parser.add_argument("--expressionProfile")
    parser.add_argument("--RSEMmodel")
    parser.add_argument("--genome", default='syn5923671')
    parser.add_argument("--gtf", default='syn7989163')
    parser.add_argument("--dipGenome", default="GRCh37v75_STAR.tar.gz")
    parser.add_argument("--inputDir", default="/home/ubuntu/inputs")
    parser.add_argument("--bucket", default="gs://smc-rna-eval")
    args = parser.parse_args()
    
    syn = synapseclient.Synapse()
    syn.login()

    expressionProfile = syn.get(args.expressionProfile)
    RSEMmodel = syn.get(args.RSEMmodel)
    genome = syn.get(args.genome)
    gtf = syn.get(args.gtf)

    dipGenome = '/'.join([args.bucket,args.dipGenome])
    subprocess.check_call(["gsutil", "cp", "-r", dipGenome, args.inputDir])
    dipGenome = '/'.join([args.inputDir,args.dipGenome])

    job = {
        "SIM_NAME" : "%s" % (args.simName),
        "GTF" : {
            "class" : "File",
            "path" : "%s" % (gtf.path)
        },
        "NUM_EVENTS" : int(args.numEvents),
        "TARGET_DEPTH" : int(args.targetDepth),
        "GENOME" : {
            "class" : "File",
            "path" : "%s" % (genome.path)
        },
        "EXPRESSION_PROFILE" : {
            "class" : "File",
            "path" : "%s" % (expressionProfile.path)
        },
        "RSEM_MODEL" : {
            "class" : "File",
            "path" : "%s" % (RSEMmodel.path)
        },
        "DIP_GENOME": {
            "class" : "File",
            "path" : "%s" % (dipGenome)
        }
    }

    input_path = os.path.join(args.inputDir, "input.json")
    with open(input_path, "w") as handle:
        handle.write(json.dumps(job, indent=4))
