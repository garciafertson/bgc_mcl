# Convert proteins from an annotated GenBank file into a FASTA file
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

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

