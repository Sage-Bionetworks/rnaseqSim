#!/usr/bin/env python

import os
from glob import glob
import synapseclient
import argparse
import yaml
import imp


def syn_login(args):
    syn = synapseclient.Synapse()
    syn.login()
    return syn

class Config:
    def __init__(self, syn, path):
        self.path = path
        with open(args.config_file) as handle:
            config = yaml.load(handle.read())
        self.config = config
        self.syn = syn
        
    def call_module(self, module, **kwargs):
        module_dir = os.path.join( os.path.dirname(os.path.abspath(self.path)), self.config[module] )

        args = {}
        for k, v in kwargs.items():
            if v.endswith("_syn"):
                args[k] = self.syn.get(self.config[v]).path
        
        print module_dir, args
        
        mod = imp.load_source('sim_runner.%s' % (module), os.path.join(module_dir, "module.py"))
        out = mod.run_module(**args)
        

def command_run(args):
    syn = syn_login(args)
    
    config = Config(syn, args.config_file)
    
    config.call_module('fusion_create', gtf_file='gtf_syn', genome_file='genome_syn')
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    subparsers = parser.add_subparsers(title="subcommand")
    
    parser_run = subparsers.add_parser('run')
    parser_run.add_argument("config_file")
    parser_run.set_defaults(func=command_run)
    
    args = parser.parse_args()
    args.func(args)

    