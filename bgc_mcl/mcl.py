'''Read porthomcl output <8.all.ort.groups> file , build jaccard index distance matrix
   and create groups using MCL algorithm
'''
import numpy as np
from collections import defaultdict
from Bio import SeqIO
import subprocess
from scripts.proteins2fasta import fasta2pfam, pfam2list
import os, subprocess

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
    a=set(a)
    b=set(b)
    ji=len(a.intersection(b))/len(a.union(b))
    return(ji)
def overlap_coef(a,b):
    a=set(a)
    b=set(b)
    oc=len(a.intersection(b))/min(len(a),len(b))
    return(oc)
def domdupsim(a,b): #Measures the amount of duplication of a shared domain
    #get elements in  union
    union=set(a).union(set(b))
    #count the number of times each domain exists in each set
    count_a=defaultdict(int)
    for element in a: #count number of times each element appears in a
        count_a[element]+=1
    count_b=defaultdict(int)
    for element in b:
        count_b[element]+=1
    #calculate S and ref Kui Lin 2006
    S=0
    for dom in union:
        S+= max(count_a[dom], count_b[dom])
    #calculate D index [0.36 ,1]
    exp=0
    D=0
    if S>0:
        for dom in union:
            x+= abs(count_a[dom]-count_b[dom])/S
        D=np.e**exp
        #D= (np.e**-x) - (x*(np.e**-x))
        #D=1-sm
    return(D)

def portho2ntwk(bgccogs, bgclist, bgcdir, output):
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
    bgc_dict= defaultdict(list)
    #dict key<-bgc_name, value<-set.add(group number)
    with open(bgclist, "r") as f:
        for infile in f:
            root_name=".".join(infile.split(".")[0:-1]) #remove .ext
            infile=root_name+".faa"
            for record in SeqIO.parse("0.input_faa/"+infile, "fasta"):
                gen_id="|".join([root_name,record.id.split("|")[2]])
                if gen_id in ortgroup_dict.keys():
                    bgc_dict[infile].append(str(ortgroup_dict[gen_id]))
                else:
                    bgc_dict[infile].append(gen_id)
    #construct jaccard index network file <network.abx>
    names=[*bgc_dict] #unpack iterable dict keys
    dim=len(names)
    network=output+".oc_network.abx"
    cutoff = 0.5 #set provisional cutoff
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
                    #print(len(a.intersection(b)))
                    out.write("%s\t%s\t%s\n" %(names[i], names[j], oc))
    #scan a range of Inflations values to run MCL algorithm on distance matrix
    return(network)

def gbk2pfam2ntwk(bgclist, outname, nohmm): #reads faa files from <0.input_faa> folder
    #open bgc file list and save list of names (remove gbk extension)
    bgcpfam_dict=defaultdict(list)
    #create "0.input_faa" folder
    try: os.makedirs(outname+"_pfam")
    except OSError: print ("Creation of the directory <output_pfam> failed")
    else: print ("Successfully created the directory <output_pfam>")
    with open(bgclist, "r") as f:
        for line in f:
            line=line.rstrip()
            root_name=".".join(line.split(".")[0:-1]) #remove .ext
            fasta_name="0.input_faa/"+root_name+".faa"
            pfam_name=outname+"_pfam/"+root_name+".dom.tbl"
            tsv_name=outname+"_pfam/"+root_name+".tsv"
            fasta2pfam(fasta_name,pfam_name)
            try: pfam2tsv(pfam_name, tsv_name) #output tsv file with pfams
            except: continue
            bgcpfam_dict[root_name]=tsv2list(tsv_name)
    names=[*bgcpfam_dict]
    dim=len(names)
    network=root_name+"pfam_network.abx"
    cutoff=0
    with open(network,"w") as out:
        for i, in range(dim):
            for j in range(i+1,dim):
                Dup_index=0
                list_a=bgcpfam_dict[names[i]]
                list_b=bgcpfam_dict[names[j]]
                if len(list_a)>0 or len(list_b)>0:
                    Dup=domdupsim(list_a,list_b)
                if Dup > cutoff:
                    out.write("%s\t%s\t%\n"%(names[i],names[j],Dup))
    return(network)

def mcl_scan(network, low, up, points, threads, output):
    for I in np.linspace(low,up,points):
        run_mcl(network, I, output, threads)

