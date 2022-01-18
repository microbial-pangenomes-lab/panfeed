#!/usr/bin/env python

import os
import sys
import logging
import pandas as pd
from pyfaidx import Fasta, Faidx
import numpy as np
import hashlib
import binascii
import shutil
from collections import namedtuple
from io import StringIO

logger = logging.getLogger('panfeed.panfeed')

def what_are_my_inputfiles(gffdir): # adds all the GFF files in the directory to a list
    filelist = []
    
    for file in os.listdir(gffdir):
        only_file = os.path.split(file)[-1]
        if file.endswith(".gff" or "gff"):
            
            genome = only_file.split(".")[0]

            logger.debug(f"Adding {genome} ({file})")

            filelist.append(genome)
        else:
            logger.debug(f"Ignoring {file}")
            
    if len(filelist) == 0:
        
        logger.error(f"No GFF files found in working directory ({gffdir})")
        
        sys.exit(1)
        
    return filelist


def prep_data_n_fasta(filelist, gffdir, output): # prepares the hybrid GFF files and creates fasta files
    
    data = {}
    
    for genome in filelist:
    
        logger.debug(f"Handling {genome}")

        logger.debug(f"Creating fasta file for {genome}")
        
        txtfile = open(os.path.join(output, f"{genome}.fasta"), "w")
    
        txtfile.write(open(os.path.join(gffdir, f"{genome}.gff"), "r").read().split("##FASTA")[1])
    
        txtfile.close()
    
        logger.debug(f"Creating faidx file for {genome}")
        
        sequences = create_faidx(os.path.join(output, f'{genome}.fasta'))
    
        logger.debug(f"Parsing features for {genome}")
        
        features = parse_gff(os.path.join(gffdir, f'{genome}.gff'))
    
        data[genome] = (sequences, features)
    
    return data

def set_input_output(stroi_in, presence_absence, output):
    # determines the names of the input, output files 
   
    logger.debug(f"Loading pangenome file from panaroo ({presence_absence})")
    genepres = pd.read_csv(presence_absence, sep=",", index_col=0).drop(columns=['Non-unique Gene name', 'Annotation'])
    # load the gene_presence_absence matrix which will be used to ascribe strain/genes to clusters
    
    if stroi_in is not None:
        logger.debug(f"Loading target strains ({stroi_in})")
        
        stroi = open(f"{stroi_in}", "r").read()
        
    else:
        logger.warning(f"No target strains provided")
        
        stroi = ""
    
    if os.path.exists(output):
       
        logger.error(f"Output directory {output} exists! Please "
                      "remove it and restart")
        sys.exit(1)
    else:
        logger.debug(f"Creating output directory ({output})")
        
        os.mkdir(output)
    
    logger.debug(f"Creating output files within {output}")
    
    kmer_stroi = open(os.path.join(output, "kmers.tsv"), "w")
    
    #creates the header for the strains of interest output file
    kmer_stroi.write("cluster\tstrain\tcontig\tcontig_start\tcontig_end\tgene_start\tgene_end\tstrand\tk-mer\n")
    
    hash_pat = open(os.path.join(output, "hashes_to_patterns.tsv"), "w")

    kmer_hash = open(os.path.join(output, "kmers_to_hashes.tsv"), "w")

    return stroi, kmer_stroi, hash_pat, kmer_hash, genepres


# a tuple with attribute names
Feature = namedtuple('Feature', ['id', 
                                 'chromosome',
                                 'start', 
                                 'end', 
                                 'strand'])

Seqinfo = namedtuple("Seqinfo", ["sequence",
                                 "chromosome",
                                 "start",
                                 "end",
                                 "strand"])
                     

Sequence = namedtuple("Sequence", [])

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
                
                logger.warning(f'{e}, skipping line "{line.rstrip()}" from {file_name}')
                
                continue
    
    return features

def iter_gene_clusters(panaroo, genome_data, up, down, patfilt):
    
    # go through each gene cluster
    for idx, row in panaroo.iterrows():
        logger.debug(f"Extracting sequences from {idx}")
        # output dict
        # key: strain
        # value: list of faidx sequence objects
        gene_sequences = {}
        
        # keep track of who has the gene
        # could be useful
        # to prepopulate the presence/absence vector
        strains = row.index
        
        sortstrain = {x: i
                        for i, x in enumerate(sorted(strains))}
        
        present = row.dropna().index
        
        absent = strains.difference(present)
        
        if patfilt == True:
            
            clusterpresab = np.zeros(len(strains)).astype(int)
            
            for strain in present:
                
                clusterpresab[sortstrain[strain]] = 1 
        else:
            
            clusterpresab = 0

        # print(clusterpresab)
        # cycle through all the strains that have the gene
        for strain, genes in row.dropna().iteritems():
            
            strain = str(strain)
            
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
                    logger.warning(f"Could not find gene {gene} from {idx}")
                    
                    continue
                # access its gene sequence
                seq = sequences[feat.chromosome][feat.start-1+up:feat.end+down]
                # be aware of strand
                revseq = -seq
                
                if feat.strand == -1:

                    revstrand = 1
                    
                else:

                    revstrand = -1
                
                # save it
                gene_sequences[strain].append(Seqinfo(str(seq), feat.chromosome, feat.start, feat.end, feat.strand))
                
                gene_sequences[strain].append(Seqinfo(str(revseq), feat.chromosome, feat.start, feat.end, revstrand))
        
        # provide an empty list if the strain does not have the gene
        for strain in absent:
            
            gene_sequences[strain] = []
            
        yield gene_sequences, idx, clusterpresab


def cluster_cutter(cluster_gen, klength, stroi, kmer_stroi, canon):
#iterates through genes present in the cluster
#shreds gene sequence into k-mers and adds them to the cluster dictionary
#if a gene belongs to a strain of interest, additional info on the k-mer is saved

    for cluster, idx, clusterpresab in cluster_gen:
        logger.debug(f"Extracting k-mers from {idx}")

        memchunk = StringIO()
        
        cluster_dict = {}
    
        sortstrain = {x: i
                      for i, x in enumerate(sorted(cluster.keys()))}
        n_strains = len(sortstrain)

        for strain in cluster.keys():

            if canon == True:#block for canonical k-mers
                
                if strain in stroi:#subblock for strains of interest
                    
                    for seq1, seq2 in zip(cluster[strain][::2], cluster[strain][1::2]):

                        num_kmer = len(seq1.sequence) - klength + 1
                        
                        seqlen = len(seq1.sequence)
                        
                        for pos in range(num_kmer):
                            
                            if pos == 0:
                            #prevents erroneous slicing of k-mers
                                revspecseq = str(seq2.sequence[-pos - klength::])
                            
                            else:
                                
                                revspecseq = str(seq2.sequence[-pos - klength:-pos:])
                            
                            contigstart = seq1.start #starting position of the gene in relation to the contig
                                
                            truestart = contigstart + pos #starting postition of the k-mer in relation to the contig
                                
                            trueend = contigstart + pos + klength #ending postition of the k-mer in relation to the contig
                            
                            specseq = str(seq1.sequence[pos:pos + klength])

                            if specseq < revspecseq: #selects the lexicographically smallest kmer
                                
                                if specseq not in cluster_dict:
                                     
                                    cluster_dict[f"{specseq}"] = np.zeros(n_strains).astype(int)
                                     
                                cluster_dict[f"{specseq}"][sortstrain[strain]] = 1

                                contig = seq1.chromosome

                                strand = seq1.strand
                                
                                if strand < 0:
                                    
                                
                                    memchunk.write(f"{idx}\t{strain}\t{contig}\t{truestart}\t{trueend}\t{pos+klength}\t{pos}\t{strand}\t{specseq}\n")
                                else:
                                    
                                    memchunk.write(f"{idx}\t{strain}\t{contig}\t{truestart}\t{trueend}\t{pos}\t{pos+klength}\t{1}\t{specseq}\n")
                                
                            else:
                                
                                strand = seq2.strand

                                contig = seq2.chromosome
                                 
                                if revspecseq not in cluster_dict:
                                     
                                    cluster_dict[f"{revspecseq}"] = np.zeros(n_strains).astype(int)
                                     
                                cluster_dict[f"{revspecseq}"][sortstrain[strain]] = 1
                                
                                contigstart = int(seq2.start)    
                                
                                if strand < 0:
                                                                    
                                    memchunk.write(f"{idx}\t{strain}\t{contig}\t{truestart}\t{trueend}\t{pos+klength}\t{pos}\t{strand}\t{revspecseq}\n")
                                
                                else:
                                                                        
                                    memchunk.write(f"{idx}\t{strain}\t{contig}\t{truestart}\t{trueend}\t{pos}\t{pos+klength}\t{1}\t{revspecseq}\n")
                                
                else:#subblock for other strains##########################################
                    
                    for seq1, seq2 in zip(cluster[strain], cluster[strain][::2]):
                        
                        num_kmer = len(seq1.sequence) - klength + 1
                        
                        seqlen = len(seq1.sequence)
                        
                        for pos, posrevcomp in zip(range(num_kmer),range(num_kmer+1,0,-1)):
                
                            specseq = str(seq1.sequence[pos:pos + klength])
                            
                            revspecseq = str(seq2.sequence[posrevcomp - klength:posrevcomp])
                            
                            canonseq = min(specseq, revspecseq)
                            
                            if canonseq not in cluster_dict:
                            
                                cluster_dict[f"{canonseq}"] = np.zeros(n_strains).astype(int)
                                     
                            cluster_dict[f"{canonseq}"][sortstrain[strain]] = 1
                
            else:#block for non-canonical k-mers#############################################################
                
                for seq in cluster[strain]:
    
                    # determines the number of k-mers in the sequence
                    num_kmer = len(seq.sequence) - klength + 1
                    
                    seqlen = len(seq.sequence)
                    
                    strand = seq.strand
                    
                    contig = seq.chromosome
                    
                    if strain in stroi:#subblock for strains of interest##########################################
                        #if a gene belongs to a strain of interest the strain,contig,start,stop,strand,sequence
                        #of its k-mers are being written to kmer_stroi
                        if strand > 0:
                        
                            for pos in range(num_kmer):
                            
                                specseq = str(seq.sequence[pos:pos + klength])
                        
                                if specseq not in cluster_dict:
                                    # init with an "empty" array
                                    cluster_dict[f"{specseq}"] = np.zeros(n_strains).astype(int)
                            
                                # flip the "0" to "1"
                                cluster_dict[f"{specseq}"][sortstrain[strain]] = 1
                                
                                contigstart = int(seq.start)
                                
                                truestart = contigstart + pos + klength
                                        
                                trueend = contigstart + pos
                        
                                memchunk.write(f"{idx}\t{strain}\t{contig}\t{truestart}\t{trueend}\t{pos}\t{pos+klength}\t{strand}\t{specseq}\n")
    
                        else:
                            
                            for pos in range(num_kmer):
                            
                                specseq = str(seq.sequence[pos:pos + klength])
                                
                                if specseq not in cluster_dict:
                                    # init with an "empty" array
                                    cluster_dict[f"{specseq}"] = np.zeros(n_strains).astype(int)
                            
                                # flip the "0" to "1"
                                cluster_dict[f"{specseq}"][sortstrain[strain]] = 1
                                
                                contigstart = int(seq.start)
                                
                                truestart = contigstart - pos
                                        
                                trueend = contigstart - pos - klength
                        
                                memchunk.write(f"{idx}\t{strain}\t{contig}\t{truestart}\t{trueend}\t{seqlen - pos}\t{seqlen - pos- klength}\t{strand}\t{specseq}\n")
                        
                    else:#subblock for other strains##########################################
                    #does the same as the if-statement block on the same indent,
                    #without logging the k-mer infos, as the strain is not a strain of interest
                        for pos in range(num_kmer):
                        
                            specseq = str(seq.sequence[pos:pos + klength])
                    
                            # specseqrev = str(seqrevcomp1[pos:pos + klength])
            
                            if specseq not in cluster_dict:
                                # init with an "empty" array
                                cluster_dict[f"{specseq}"] = np.zeros(n_strains).astype(int)
                        
                            # flip the "0" to "1"
                            cluster_dict[f"{specseq}"][sortstrain[strain]] = 1
                        
                        
        kmer_stroi.write(memchunk.getvalue())
            
        yield cluster_dict, clusterpresab

        

def pattern_hasher(cluster_dict_iter, hash_pat, kmer_hash, genepres, patfilt):
    #iterates through the cluster dictionary output by cluster_cutter()
    #outputs the hashed patterns, patterns and k-mers to the output files
    #two files are created: hashed k-mer patterns to presence/absence patterns
    #                       k-mer sequences to hashed k-mer patterns
    memchonkheader1 = StringIO()
    
    memchonkheader2 = StringIO()
    
    memchonkheader1.write("hashed_pattern")

    for strain in sorted(genepres.columns):
        
        memchonkheader1.write(f"\t{strain}")
    
    memchonkheader1.write("\n")
    
    hash_pat.write(memchonkheader1.getvalue())
    
    memchonkheader2.write("k-mer\thashed_pattern\n")
    
    kmer_hash.write(memchonkheader2.getvalue())

    for cluster_dict, clusterpresab in cluster_dict_iter:
        
        memchunkhash_pat = StringIO()
        
        memchunkkmer_hash = StringIO()
        
        if patfilt == True:
            
            for kmer in cluster_dict:
    
                if tuple(cluster_dict[kmer]) != tuple(clusterpresab):
                    
                    pattern = cluster_dict[kmer].view(np.uint8)
    
                    hashed = hashlib.md5(pattern)     
                    
                    khash = binascii.b2a_base64(hashed.digest()).decode()[:24]
                                
                    patterntup = "\t".join(map(str, cluster_dict[kmer]))
                                    
                    memchunkhash_pat.write(f"{khash}\t{patterntup}\n")
                                    
                    memchunkkmer_hash.write(f"{kmer}\t{khash}\n")
                   
        else:
                    
            for kmer in cluster_dict:
                    
                pattern = cluster_dict[kmer].view(np.uint8)
                        
                hashed = hashlib.md5(pattern)
                            
                khash = binascii.b2a_base64(hashed.digest()).decode()[:24]
                        
                patterntup = "\t".join(map(str, cluster_dict[kmer]))
                            
                memchunkhash_pat.write(f"{khash}\t{patterntup}\n")
                            
                memchunkkmer_hash.write(f"{kmer}\t{khash}\n")
            
        hash_pat.write(memchunkhash_pat.getvalue())
        
        kmer_hash.write(memchunkkmer_hash.getvalue())


if __name__ == "__main__":
    pass
