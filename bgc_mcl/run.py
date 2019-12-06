#!/usr/bin/env python
import subprocess
import os, sys
from bgc_mcl.porthomcl import porthomcl_analysis
from bgc_mcl.mcl import mcl_scan
class Run:
    def __init__(self, args):
        self.args=args
    def porthomcl(self): 
        bgclist=self.args.list
        bgcdir=self.args.dir
        threads=self.args.threads
        run_porthomcl(bgclist, bgcdir, threads)
    def mcl(self):
        bgclist=self.args.list
        bgcdir=self.args.dir
        bgccogs=self.args.bgccogs
        low=self.args.inf_low
        up=self.args.inf_upper
        inc=self.args.increment
        outname=self.args.output
        mcl_scan(bgccogs, bgclist, bgcdir, low, up, inc, outname)
        
    def main(self):
        if self.args.subparser_name=="porthomcl":
            self.porthomcl()
            print("running POrthoMCL on BGC.gbk files\n")
        elif self.args.subparser_name=="mcl":
            self.mcl()
            print("running MCL on BGC.gbk using %s BGCcogs" %self.args.)



