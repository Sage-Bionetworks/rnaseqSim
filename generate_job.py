#!/usr/bin/env python

import os
import yaml
import argparse
import json
import synapseclient

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--simName")
    parser.add_argument("--numEvents")
    parser.add_argument("--targetDepth")
    parser.add_argument("--expressionProfile")
    parser.add_argument("--RSEMmodel")
    #parser.add_argument("--genome")
    #parser.add_argument("--gtf")    

    args = parser.parse_args()
    
    syn = synapseclient.Synapse()
    syn.login()

    expressionProfile = syn.get(args.expressionProfile)
    RSEMmodel = syn.get(args.RSEMmodel)
    genome = syn.get('syn5923671')
    gtf = syn.get('syn7989163')

    job = {
        "SIM_NAME" : "%s" % (args.simName),
        "GTF" : {
            "class" : "File",
            "path" : "%s" % (gtf.path)
        },
        "NUM_EVENTS" : args.numEvents,
        "TARGET_DEPTH" : args.targetDepth,
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
        }
    }
        
    print json.dumps(job)
