import os, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

# Convert proteins from an annotated GenBank file into a FASTA file
def gbk2faa(infile, outfile):
    with open (outfile, 'w') as out:
        print("converting gbk into faa", infile)
        for record in SeqIO.parse(infile, "genbank"):
            print('Converting proteins from {}...'.format(record.id))
            proteins = []
            for feature in record.features:
                if feature.type == 'CDS' and feature.qualifiers.get('translation'):
                    locus_tags = feature.qualifiers.get('locus_tag')
                    locus_tag = str(locus_tags[0]) if locus_tags else 'unknown_locus_tag'
                    # a unique protein id is required, if protein_id is not present, we need a locus tag
                    location = str(feature.location.start) + '-' + str(feature.location.end)
                    protein_ids = feature.qualifiers.get('protein_id')
                    if protein_ids:
                        protein_id = str(protein_ids[0])
                    elif locus_tags:
                        protein_id = str(locus_tags[0])
                    else:
                        protein_id = 'unknown_protein'+location
                    r_id = record.id + '|' + locus_tag + '|' + protein_id + '|' + location + '|' + str(feature.location.strand)
                    translation = Seq(feature.qualifiers.get('translation')[0], generic_protein)
                    r = SeqRecord(translation, r_id, description='')
                    proteins.append(r)
        SeqIO.write(proteins, out, 'fasta')
    print('Saved faa file', outfile)
# run hmmscan into faa file, save PFAM env varialble with path tho Pfam database

def fasta2pfam(infile, outfile):
    #Run hmmscan using fasta file input output names domain table 
    try:
        path_pfamdb=str(os.environ['PFAM_A']) #enviroment variable with path PFAM database
    except OSError:
        print("Please set enviromental varialble PFAM_A with path to pfam database")
    cmd=["hmmscan", "--domtblout", outfile, "--cut_tc","--noali", path_pfamdb, infile]
    try:
        subprocess.call(cmd)
    except OSError:
        print("mcl command failed, check whether enviroment variable PFAM_A is defined")

# convert hmmscan domain table output into pfam tsv and pfam list for bgc file.
from Bio import SearchIO
import pandas as pd
import numpy as np
def pfam2list(infile):
    """
    Convert hmmscan tabular format into python list
    :infile: Path to HMMscan tabular format result file
    :return: domain list and print Domain DataFrame.
    """
    # Read domain matches in all proteins
    queries = SearchIO.parse(infile, 'hmmscan3-domtab')
    # Extract all matched domain hits
    print('Reading domains...')
    domains = []
    for query in queries:
        query_domains = []
        for hit in query.hits:
            best_index = np.argmin([hsp.evalue for hsp in hit.hsps])
            best_hsp = hit.hsps[best_index]
            pfam_id = hit.accession.split('.')[0]
            query_domains.append({
                'pfam_id': pfam_id,
                'sequence_id': query.id,
                'domain_start': int(best_hsp.hit_start),
                'domain_end': int(best_hsp.hit_end),
                'evalue': float(best_hsp.evalue),
                'bitscore' : float(best_hsp.bitscore)
            })
        domains += sorted(query_domains, key=lambda x: x['domain_start'])

    print('Domain hits total: {}'.format(len(domains)))
    return(domains)

