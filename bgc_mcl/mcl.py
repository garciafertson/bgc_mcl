'''Read porthomcl output <8.all.ort.groups> file , build jaccard index distance matrix
   and create groups using MCL algorithm
'''
import numpy as np
from collections import defaultdict
from Bio import SeqIO
import subprocess

def run_mcl(infile, inflation, outfile,threads):
    out=outfile+".I_"+str(inflation)[0:3]+".txt"
    cmd=["mcl", infile, "-te", str(threads), "-I", str(inflation), "-o", out, "--abc"]
    try:
        subprocess.call(cmd)
    except OSError:
        print("mcl command failed, chech if mcl program is installed")
    else:
        print("Running mcl program, Inflation value:", inflation)

def jaccard_index(a,b):
    #calculate jaccard index for pairs of bgc
    ji=len(a.intersection(b))/len(a.union(b))
    return(ji)
def overlap_coef(a,b):
    oc=len(a.intersection(b))/min(len(a),len(b))
    return(oc)

def mcl_scan(bgccogs, bgclist, bgcdir, low, up, points, output,threads):
    #read groups of othologues from bgccogs 8.all.oth.group file
    #and create dict key<-fasta_id, value <-line_number or group number
    ortgroup_dict={}
    i=0
    with open (bgccogs, "r") as f:
        for line in f:
            line=line.rstrip()
            for fasta_id in line.split("\t"):
                ortgroup_dict[fasta_id]=i
            i+=1
    #open list of gbk files, read corresponding fasta files
    #retrieve id and define set of BGCcog for each fasta file
    bgc_dict= defaultdict(set)
    #dict key<-bgc_name, value<-set.add(group number)
    with open(bgclist, "r") as f:
        for infile in f:
            root_name=infile.split(".")[0]
            infile=root_name+".faa"
            for record in SeqIO.parse("0.input_faa/"+infile, "fasta"):
                gen_id="|".join([root_name,record.id.split("|")[2]])
                if gen_id in ortgroup_dict.keys():
                    bgc_dict[infile].add(str(ortgroup_dict[gen_id]))
                else:
                    bgc_dict[infile].add(gen_id)
    #construct jaccard index network file <network.abx>
    names=[*bgc_dict] #unpack iterable dict
    dim=len(names)
    network=output+".oc_network.abx"
    cutoff = 0 #set provisional cutoff
    with open(network,"w") as out:
        for i in range(dim):
            for j in range(i+1,dim):
                oc=ji=0
                a=bgc_dict[names[i]]
                b=bgc_dict[names[j]]
                if len(a)>0 and len(b)>0:
                    ji=jaccard_index(a,b)
                    oc=overlap_coef(a,b)
                if oc > cutoff:
                    print(len(a.intersection(b)))
                    out.write("%s\t%s\t%s\n" %(names[i], names[j], oc))
    #scan a range of Inflations values to run MCL algorithm on distance matrix
    for I in np.linspace(low,up,points):
        run_mcl(network, I, output, threads)
