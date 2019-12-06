'''
Recieve filelist of gbk BGC and runs PorthoMCL scaning a range of different inflation values 
'''
import os, subprocess
from scripts.proteins2fasta import gbk2faa

# run porthoMCL
def run_porthomcl(threads): #folder contains the path to BGC.gbk
    path_porhtomcl=str(os.environ['PORTHOMCL']) enviroment variable with path to external program
    cmd=[path_porthomcl+"porthomcl.sh", "-t", threads]
    try:subprocess.call(cmd)
    except: print("Plase set enviromental variable PORTHOMCL with path to <porthomcl.sh>")
    else:print("Running porthomcl.sh, please be patient")

def porthomcl_analysis(bgclist, bgcdir, threads):
    #create "0.input_faa" folder
    try: os.makedirs("0.input_faa")
    except OSError: print ("Creation of the directory 0.input_faa failed")
    else: print ("Successfully created the directory 0.input_faa")
    #convert genbank files into fasta files 
    with open(bgclist, "r") as f:
        for line in f:
            root=line.split(".")[0]
            in_gbk=bgcdir+line
            out_faa="0.input_faa/"+root+".faa"
            gbk2faa(in_gbk, out_faa)
    #run portomcl.sh
    run_porthomcl(threads)

