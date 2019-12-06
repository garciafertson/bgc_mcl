'''Read porthomcl output <8.all.ort.groups> file , build jaccard index distance matrix
   and create groups using MCL algorithm
'''
import numpy as np
import defaultdict
from Bio import seqIO

def run_mcl(matrix, inflation):
    cmd=["mcl", matrix, "-te", threads, "-I", inflation]
    try:subprocess.call(cmd)
    except: print("mcl command failed, chech if mcl program is installed")
    else:print("Running mcl program, Inflation value:", inflation)

def jaccard_index(a,b):
    #calculate jaccard index for pairs of bgc
    ji=len(a.intersection(b))/len(a.union(b))
    return(ji)
def overlap_coef(a,b):
    oc=len(a.intersection(b))/min(len(a),len(b))
    return(oc)

def mcl_scan(bgccogs, bgclist, bgcdir, low, up, inc, output):
    #read groups of othologues from bgccogs 8.all.oth.group file
    #and create dict key<-fasta_id, value <-line_number or group number
    ortgroup_dict={}
    i=0
    with open (bgccogs, "r") as f:
        for line in f:
            for fasta_id in line.split("\t"):
                ortgroup_dic[fasta_id]=i
            i+=1
    #read fasta files, retrieve id
    bgc_dict= defaultdict(set)
    #and create dict key<-bgc_name, value<-set.add(group number)
    with open(bgclist, "r") as f:
        for line in f:


    #construct jaccard index distance matrix
    matrix=np.zeros(len())
    with open(bgclist,"r") as f:
        for line in f:
            bgcname=line.split(".")[0]

    #output <graph.mci> or <input.abx>

    #scan a range of Inflations values to run MCL algorithm on distance matrix
    for i in range(low, up, inc):
        run_mcl(out_matrix, i)
