#!/usr/bin/env python
import subprocess
import os, sys
from bgc_mcl.porthomcl import porthomcl_analysis
from bgc_mcl.mcl import portho2ntwk, gbk2pfam2ntwk, mcl_scan 


class Run:
    def __init__(self, args):
        self.args=args
    def porthomcl(self): 
        bgclist=self.args.list
        bgcdir=self.args.dir
        threads=self.args.threads
        porthomcl_analysis(bgclist, bgcdir, threads)
    def mcl(self):
        bgclist=self.args.list
        bgcdir=self.args.dir
        bgccogs=self.args.bgccogs
        low=self.args.inf_lower
        up=self.args.inf_upper
        points=self.args.points
        outname=self.args.output
        threads=self.args.threads
        cutoff=self.args.cutoff
        index=self.args.index
        if cutoff < 0 and cutoff > 1: sys.exit("Error: --cutoff value out of range [0,1]")
        network=portho2ntwk(bgccogs, bgclist, bgcdir, outname, sim_index, cutoff)
        mcl_scan(network,low,up,points,threads,outname)
    def pfammcl(self):
        bgclist=self.args.list
        bgcdir=self.args.dir
        outname=self.args.output
        low=self.args.inf_lower
        up=self.args.inf_upper
        points=self.args.points
        threads=self.args.threads
        network=gbk2pfam2ntwk(bgclist, outname)
        mcl_scan(network, low,up,points, threads)

    def main(self):
        if self.args.subparser_name=="porthomcl":
            self.porthomcl()
            print("running Parallel OrthoMCL on BGC.gbk files\n")
        elif self.args.subparser_name=="mcl":
            self.mcl()
            print("running MCL on BGC.gbk using %s BGCcogs" %self.args.list)
        elif self.args.subparser_name=="pfammcl":
            self.pfammcl()
            print("runing MCL clustering on BGC.gbk files %s using PFAM content" %self.args.list)


