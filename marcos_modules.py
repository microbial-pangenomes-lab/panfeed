import os
import sys
import pandas as pd
from pyfaidx import Fasta
from collections import namedtuple
from io import StringIO
# a tuple with attribute names
Feature = namedtuple('Feature', ['id', 
                                 'chromosome',
                                 'start', 
                                 'end', 
                                 'strand'])

def create_faidx(file_name, uppercase=True):
    faidx = Fasta(file_name,
                  sequence_always_upper=uppercase)
    return faidx

def parse_gff(file_name, feature_types=None):
    if feature_types is None:
        feature_types = {'CDS'}

    # output dict
    # key: feature ID
    # value: Feature NamedTuple
    features = {}
    with open(file_name, 'r') as gff:
        for line in gff:
            if line.lstrip().startswith('##FASTA'):
                # start of FASTA entries, end of file
                break
            elif line.lstrip().startswith('#'):
                # comment, ignore
                continue
            
            # should be a valid GFF3 line
            entries = line.split('\t')
            
            try:
                ftype = entries[2]
                if ftype not in feature_types:
                    continue

                chrom = entries[0]
                start = int(entries[3])
                end = int(entries[4])
                strand = entries[6]
                # integer takes up less space
                if strand == '+':
                    strand = 1
                else:
                    strand = -1

                # fetch the feature ID from the last field
                ID = None
                for entry in entries[8].split(';'):
                    if entry.startswith('ID') and '=' in entry:
                        ID = entry.split('=')[1]
                # could not find it, skip this entry
                if ID is None:
                    continue

                # save the relevant details
                features[ID] = Feature(ID, chrom, start, end, strand)
            except Exception as e:
                # not distinguishing between exceptions
                # not great behaviour
                sys.stderr.write(f'{e}, skipping line "{line.rstrip()}"\n')
                continue
    
    return features

def iter_gene_clusters(panaroo, genome_data, error_log, up, down):
    errormem = StringIO()
    # go through each gene cluster
    for idx, row in panaroo.iterrows():
        # output dict
        # key: strain
        # value: list of faidx sequence objects
        gene_sequences = {}
        
        # keep track of who has the gene
        # could be useful
        # to prepopulate the presence/absence vector
        strains = row.index
        present = row.dropna().index
        absent = strains.difference(present)

        # cycle through all the strains that have the gene
        for strain, genes in row.dropna().iteritems():
            gene_sequences[strain] = []
            # access the GFF/fasta data
            sequences, features = genome_data[strain]
            # be aware of paralogs
            for gene in genes.split(';'):
                # access a particular gene
                try:
                    feat = features[gene]
                except KeyError:
                    # e.g. refound genes are not in the GFF
                    errormem.write(f'Could not find gene {gene}\n')
                    continue
                # access its gene sequence
                seq = sequences[feat.chromosome][feat.start-1+up:feat.end+down]
                # be aware of strand
                if feat.strand < 0:
                    seq = -seq
                    
                # save it
                gene_sequences[strain].append(seq)
        
        # provide an empty list if the strain does not have the gene
        for strain in absent:
            gene_sequences[strain] = []
        
        yield gene_sequences, idx
        
    error_log.write(errormem.getvalue())
    
# data = {}

# # hard coded for this example
# genomes = ['12001', '12020', '12022']
# for genome in genomes:
#     sequences = create_faidx(f'{genome}.fasta')
#     features = parse_gff(f'{genome}.gff')
    
#     data[genome] = (sequences, features)

# # read panaroo's output
# pg = pd.read_csv('gene_presence_absence.csv', sep=',',
#                  index_col=0)
# # drop useless columns
# pg = pg.drop(columns=['Non-unique Gene name', 'Annotation'])

# # just go through the first gene
# for gene in iter_gene_clusters(pg, data):
#     break

# for each gene we get this dictionary
# key: strain
# value: faidx sequence object (contains location info also)


# {'12001': [>NC_000913.3:1531816-1532952
  # ATGGAACTTAAAAAATTGATGGGACATATTTCTATTATCCCCGATTACAGACAAGCCTGGAAAATGGAACATAAGTTATCGGATATTCTACTGTTGACTATTTGTGCCGTTATTTCTGGTGCAGAAGGCTGGGAAGATATAGAGGATTTTGGGGAAACACATCCCGATTTTTTGAAGCAATATGGTGATTTTGAAAATGGTATTCCTGTTCACGACACCATTGCCAGAGTTGTATCCTGTATCAGTCCTGCAAAATTTCACGAGTGCTTTATTAACTGGATGCGTGACTGCCATTCTTCAGATGATAAAGACGTCATTGCAATTGATGGAAAAACGCTCCGGCATTCTTATGATAAGAGTCGCCGCAGGGGAGCGATTCATGTCATTAGTGCGTTCTCAACAATGCACAGTCTGGTCATCGGACAGATCAAGACGGATGAGAAATCTAATGAGATTACAGCTATCCCAGAACTTCTTAACATGCTGGATATTAAAGGAAAAATCATCACAACTGATGCGATGGGTTGCCAGAAAGATATTGCAGAGAAGATACAAAAACAGGGAGGTGATTATTTATTCGCGGTAAAAGGAAACCAGGGGCGGCTAAATAAAGCCTTTGAGGAAAAATTTCCGCTGAAAGAATTAAATAATCCAGCGCATGACAGTTACGCAATGAGTGAAAAGAGTCACGGCAGAGAAGAAATCCGTCTTCATATTGTTTGCGATGTCCCTGATGAACTTATTGATTTCACGTTTGAATGGAAAGGGCTGAAGAAATTATGCGTGGCAGTCTCCTTTCGGTCCATAATAGCAGAACAAAAGAAAGAGCTCGAAATGACGGTCAGATATTATATCAGTTCTGCTGATTTAACCGCAGAGAAGTTCGCCACAGCAATCCGAAACCACTGGCATGTGGAGAATAAGCTGCACTGGCGTCTGGACGTGGTAATGAATGAAGACGACTGCAAAATAAGAAGAGGAAATGCAGCAGAATTATTTTCAGGGATACGGCACATTGCTATTAATATTTTGACGAATGATAAGGTATTCAAGGCAGGGTTAAGACGTAAGATGCGAAAAGCAGCCATGGACAGAAACTACCTGGCGTCAGTCCTTACGGGGAGCGGGCTTTCGTAA],
  # '12020': [>ECOR-04_contig000003:45-347
  # ATGACGGTCAGATATTATATCAGTTCTGCTGATTTAACCGCAGAAAAGTTCGCCACAGCAATCCGAAACCACTGGCACGTGGAGAATAAGCTGCACTGGCGTATGGACGTGGTAATGAATGAAGACGACTGCAAAATAAGAAGAGGAAACGCCGCAGAATTATTTTCAGGGATACGGCACATCGCTATTAATATTTTAACGAATGATAAGGTATTCAAGGCAGGGTTAAGACGTAAGATGCGAAAAGCAGCCATGGATAGAAACTATCTCGCGTCAGTCCTTGCGGGGAGCGGGCTTTCGTAA],
  # '12022': [>ECOR-06_contig000023:83-385
  # ATGACGGTCAGATATTATATCAGTTCTGCTGATTTAACCGCTGAGAAATTCGCCACAGCGATCCGAAATCACTGGCACGTGGAGAATAAGCTGCACTGGCGTCTGGACGTGGTAATGAATGAAGACGACTGCAAAATAAGAAGAGGAAATGCAGCAGAATTATTTTCAGGGATACGGCACATTGCTATTAATATTTTGACGAATGATAAGGTATTCAAGGCAGGGTTAAGACGTAAGATGCGAAAAGCAGCCATGGACAGAAACTATCTGGCGTCAGTCCTTGCGGGGAGCGGGCTTTCGTAA]}

# pick an example faidx sequence object
# look at its attributes
# seq = gene['12001'][0]

# # get k-mers out of me!
# seq.seq

# 'ATGGAACTTAAAAAATTGATGGGACATATTTCTATTATCCCCGATTACAGACAAGCCTGGAAAATGGAACATAAGTTATCGGATATTCTACTGTTGACTATTTGTGCCGTTATTTCTGGTGCAGAAGGCTGGGAAGATATAGAGGATTTTGGGGAAACACATCCCGATTTTTTGAAGCAATATGGTGATTTTGAAAATGGTATTCCTGTTCACGACACCATTGCCAGAGTTGTATCCTGTATCAGTCCTGCAAAATTTCACGAGTGCTTTATTAACTGGATGCGTGACTGCCATTCTTCAGATGATAAAGACGTCATTGCAATTGATGGAAAAACGCTCCGGCATTCTTATGATAAGAGTCGCCGCAGGGGAGCGATTCATGTCATTAGTGCGTTCTCAACAATGCACAGTCTGGTCATCGGACAGATCAAGACGGATGAGAAATCTAATGAGATTACAGCTATCCCAGAACTTCTTAACATGCTGGATATTAAAGGAAAAATCATCACAACTGATGCGATGGGTTGCCAGAAAGATATTGCAGAGAAGATACAAAAACAGGGAGGTGATTATTTATTCGCGGTAAAAGGAAACCAGGGGCGGCTAAATAAAGCCTTTGAGGAAAAATTTCCGCTGAAAGAATTAAATAATCCAGCGCATGACAGTTACGCAATGAGTGAAAAGAGTCACGGCAGAGAAGAAATCCGTCTTCATATTGTTTGCGATGTCCCTGATGAACTTATTGATTTCACGTTTGAATGGAAAGGGCTGAAGAAATTATGCGTGGCAGTCTCCTTTCGGTCCATAATAGCAGAACAAAAGAAAGAGCTCGAAATGACGGTCAGATATTATATCAGTTCTGCTGATTTAACCGCAGAGAAGTTCGCCACAGCAATCCGAAACCACTGGCATGTGGAGAATAAGCTGCACTGGCGTCTGGACGTGGTAATGAATGAAGACGACTGCAAAATAAGAAGAGGAAATGCAGCAGAATTATTTTCAGGGATACGGCACATTGCTATTAATATTTTGACGAATGATAAGGTATTCAAGGCAGGGTTAAGACGTAAGATGCGAAAAGCAGCCATGGACAGAAACTACCTGGCGTCAGTCCTTACGGGGAGCGGGCTTTCGTAA'

# # also me!
# seq.complement.seq

# 'TACCTTGAATTTTTTAACTACCCTGTATAAAGATAATAGGGGCTAATGTCTGTTCGGACCTTTTACCTTGTATTCAATAGCCTATAAGATGACAACTGATAAACACGGCAATAAAGACCACGTCTTCCGACCCTTCTATATCTCCTAAAACCCCTTTGTGTAGGGCTAAAAAACTTCGTTATACCACTAAAACTTTTACCATAAGGACAAGTGCTGTGGTAACGGTCTCAACATAGGACATAGTCAGGACGTTTTAAAGTGCTCACGAAATAATTGACCTACGCACTGACGGTAAGAAGTCTACTATTTCTGCAGTAACGTTAACTACCTTTTTGCGAGGCCGTAAGAATACTATTCTCAGCGGCGTCCCCTCGCTAAGTACAGTAATCACGCAAGAGTTGTTACGTGTCAGACCAGTAGCCTGTCTAGTTCTGCCTACTCTTTAGATTACTCTAATGTCGATAGGGTCTTGAAGAATTGTACGACCTATAATTTCCTTTTTAGTAGTGTTGACTACGCTACCCAACGGTCTTTCTATAACGTCTCTTCTATGTTTTTGTCCCTCCACTAATAAATAAGCGCCATTTTCCTTTGGTCCCCGCCGATTTATTTCGGAAACTCCTTTTTAAAGGCGACTTTCTTAATTTATTAGGTCGCGTACTGTCAATGCGTTACTCACTTTTCTCAGTGCCGTCTCTTCTTTAGGCAGAAGTATAACAAACGCTACAGGGACTACTTGAATAACTAAAGTGCAAACTTACCTTTCCCGACTTCTTTAATACGCACCGTCAGAGGAAAGCCAGGTATTATCGTCTTGTTTTCTTTCTCGAGCTTTACTGCCAGTCTATAATATAGTCAAGACGACTAAATTGGCGTCTCTTCAAGCGGTGTCGTTAGGCTTTGGTGACCGTACACCTCTTATTCGACGTGACCGCAGACCTGCACCATTACTTACTTCTGCTGACGTTTTATTCTTCTCCTTTACGTCGTCTTAATAAAAGTCCCTATGCCGTGTAACGATAATTATAAAACTGCTTACTATTCCATAAGTTCCGTCCCAATTCTGCATTCTACGCTTTTCGTCGGTACCTGTCTTTGATGGACCGCAGTCAGGAATGCCCCTCGCCCGAAAGCATT'

# # chromosome, start, end, strand
# seq.name, seq.start, seq.end, seq.orientation

# ('NC_000913.3', 1531816, 1532952, 1)

