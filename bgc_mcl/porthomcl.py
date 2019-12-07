'''
Recieve filelist of gbk BGC and runs PorthoMCL scaning a range of different inflation values 
'''
import os, subprocess
from scripts.proteins2fasta import gbk2faa

# run porthoMCL
def run_porthomcl(threads): #folder contains the path to BGC.gbk
    try:
        path_porthomcl=str(os.environ['PORTHOMCL']) #enviroment variable with path to external program
    except OSError: 
        print("Please set enviromental varialble PORTHOMCL=/path/to/porthomcl/")
    program=path_porthomcl+"/porthomcl.sh"
    print(program)
    threads=str(threads)
    cmd=["bash", program, "-t", threads, "-f", "3"]
    print("Running PorthoMCL, please wait...")
    try:
        subprocess.call(cmd)
    except OSError:
        print("Program <porthomcl.sh> not found")

def porthomcl_analysis(bgclist, bgcdir, threads):
    #create "0.input_faa" folder
    try: os.makedirs("0.input_faa")
    except OSError: print ("Creation of the directory 0.input_faa failed")
    else: print ("Successfully created the directory 0.input_faa")
    #convert genbank files into fasta files 
    with open(bgclist, "r") as f:
        for line in f:
            line=line.rstrip()
            root=line.split(".")[0]
            in_gbk=bgcdir+line
            out_faa="0.input_faa/"+root+".faa"
            gbk2faa(in_gbk, out_faa)
    #run portomcl.sh
    run_porthomcl(threads)

