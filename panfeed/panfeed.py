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


def clean_up_fasta(filelist, output):
    for genome in filelist:
        fasta_file = os.path.join(output, f"{genome}.fasta")
        logger.debug(f"Removing fasta file for {genome} ({fasta_file})")
        if os.path.isfile(fasta_file):
            os.remove(fasta_file)
        else:
            logger.warning(f"Could not delete {fasta_file}")
        
        faidx_file = os.path.join(output, f"{genome}.fasta.fai")
        logger.debug(f"Removing faidx file for {genome} ({faidx_file})")
        if os.path.isfile(faidx_file):
            os.remove(faidx_file)
        else:
            logger.warning(f"Could not delete {faidx_file}")

def set_input_output(stroi_in, presence_absence, output, single_file=True):
    # determines the names of the input, output files 
   
    logger.debug(f"Loading pangenome file from panaroo ({presence_absence})")
    genepres = pd.read_csv(presence_absence,
                           sep=",", index_col=0,
                           low_memory=False).drop(
                                   columns=['Non-unique Gene name', 'Annotation'])
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
    
    if single_file:
        kmer_stroi = create_kmer_stroi(output)
        hash_pat, kmer_hash = create_hash_files(output) 
    else:
        # delayied opening of the files
        kmer_stroi = None
        hash_pat = None
        kmer_hash = None

    return stroi, kmer_stroi, hash_pat, kmer_hash, genepres

def create_kmer_stroi(output):
    kmer_stroi = open(os.path.join(output, "kmers.tsv"), "w")
    
    #creates the header for the strains of interest output file
    kmer_stroi.write("cluster\tstrain\tfeature_id\tcontig\tcontig_start\tcontig_end\tgene_start\tgene_end\tstrand\tk-mer\n")

    return kmer_stroi

def create_hash_files(output):
    hash_pat = open(os.path.join(output, "hashes_to_patterns.tsv"), "w")
    kmer_hash = open(os.path.join(output, "kmers_to_hashes.tsv"), "w")

    return hash_pat, kmer_hash


# a tuple with attribute names
Feature = namedtuple('Feature', ['id', 
                                 'chromosome',
                                 'start', 
                                 'end', 
                                 'strand'])

Seqinfo = namedtuple("Seqinfo", ["sequence",
                                 "compsequence",
                                 "id",
                                 "chromosome",
                                 "start",
                                 "end",
                                 "strand",
                                 "offset"])
                     

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

def iter_gene_clusters(panaroo, genome_data, up, down, down_start_codon, patfilt):
    
    # go through each gene cluster
    all_ogs = panaroo.shape[0]
    for i, (idx, row) in enumerate(panaroo.iterrows()):
        logger.debug(f"Extracting sequences from {idx} ({i+1}/{all_ogs})")
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
                    logger.warning(f"Could not find gene {gene} from {idx} in {strain}")
                    
                    continue
                
                # corner case: upstream offset is over the contig's edge
                # adjust if so and save the true offset
                if feat.strand > 0 and feat.start-1-up < 0:
                    offset = feat.start - 1
                elif feat.strand < 0 and feat.end+up > len(sequences[feat.chromosome]):
                    offset = len(sequences[feat.chromosome]) - feat.end
                else:
                    offset = up
                
                # access its gene sequence
                if not down_start_codon:
                    # down means relative to the gene's 3'
                    if feat.strand > 0:
                        seq = sequences[feat.chromosome][feat.start-1-offset:feat.end+down]
                        seq_start = feat.start - offset
                        seq_end = feat.end + down
                    else:
                        seq = -sequences[feat.chromosome][feat.start-1-down:feat.end+offset]
                        seq_start = feat.start - down
                        seq_end = feat.end + offset
                else:
                    # down means relative to the gene's 5'
                    if feat.strand > 0:
                        seq = sequences[feat.chromosome][feat.start-1-offset:feat.start+down]
                        seq_start = feat.start - offset
                        seq_end = feat.start + down
                    else:
                        seq = -sequences[feat.chromosome][feat.end-1-down:feat.end+offset]
                        seq_start = feat.end - down
                        seq_end = feat.end + offset

                # reverse complement
                revseq = -seq
                # reverse it, easier to slice
                # so effectively now it's just the complement
                revseq = revseq[::-1]

                # save it
                seq1 = Seqinfo(str(seq), str(revseq),
                               feat.id, feat.chromosome,
                               seq_start, seq_end,
                               feat.strand, offset)
                gene_sequences[strain].append(seq1)
        
        # provide an empty list if the strain does not have the gene
        for strain in absent:
            
            gene_sequences[strain] = []
            
        yield gene_sequences, idx, clusterpresab


def cluster_cutter(cluster_gen, klength, stroi, kmer_stroi, canon, output):
#iterates through genes present in the cluster
#shreds gene sequence into k-mers and adds them to the cluster dictionary
#if a gene belongs to a strain of interest, additional info on the k-mer is saved

    # set flag for multiple files
    multiple_files = False
    if kmer_stroi is None:
        multiple_files = True

    for cluster, idx, clusterpresab in cluster_gen:
        logger.debug(f"Extracting k-mers from {idx}")

        memchunk = StringIO()
        if multiple_files:
            # close existing file?
            if kmer_stroi is not None:
                kmer_stroi.close()
            # we need to create the file
            path = os.path.join(output, idx)
            if not os.path.exists(path):
                os.mkdir(path)
            kmer_stroi = create_kmer_stroi(path)
        
        cluster_dict = {}
    
        sortstrain = {x: i
                      for i, x in enumerate(sorted(cluster.keys()))}
        n_strains = len(sortstrain)

        for strain in cluster.keys():
            for seq in cluster[strain]:
                gene_id = seq.id
                offset = seq.offset

                num_kmer = len(seq.sequence) - klength + 1
                seqlen = len(seq.sequence)

                contig = seq.chromosome
                strand = seq.strand
                for pos in range(num_kmer):
                    specseq = str(seq.sequence[pos:pos + klength])
                    # get back to reverse complement
                    revspecseq = str(seq.compsequence[pos:pos + klength])[::-1]
                
                    if canon == True:
                        if specseq <= revspecseq:
                            canonseq = specseq
                            used_strand = strand
                        else:
                            canonseq = revspecseq
                            used_strand = - strand

                        if canonseq not in cluster_dict:
                            cluster_dict[f"{canonseq}"] = np.zeros(n_strains).astype(int)
                        cluster_dict[f"{canonseq}"][sortstrain[strain]] = 1
                    else:
                        used_strand = strand
                        if specseq not in cluster_dict:
                            cluster_dict[f"{specseq}"] = np.zeros(n_strains).astype(int)
                        cluster_dict[f"{specseq}"][sortstrain[strain]] = 1
                        
                        if revspecseq not in cluster_dict:
                            cluster_dict[f"{revspecseq}"] = np.zeros(n_strains).astype(int)
                        cluster_dict[f"{revspecseq}"][sortstrain[strain]] = 1

                    if strain in stroi:
                        if strand > 0:
                            contigstart = seq.start
                            truestart = contigstart + pos
                            trueend = contigstart + pos + klength

                        else:
                            contigstart = seq.end
                            trueend = contigstart - pos
                            truestart = contigstart - pos - klength

                        genestart = pos - offset
                        geneend = pos + klength - offset
                        if canon == True:
                            memchunk.write(f"{idx}\t{strain}\t{gene_id}\t{contig}\t{truestart}\t{trueend}\t{genestart}\t{geneend}\t{used_strand}\t{canonseq}\n")
                        else:
                            memchunk.write(f"{idx}\t{strain}\t{gene_id}\t{contig}\t{truestart}\t{trueend}\t{genestart}\t{geneend}\t{used_strand}\t{specseq}\n")
                            memchunk.write(f"{idx}\t{strain}\t{gene_id}\t{contig}\t{truestart}\t{trueend}\t{genestart}\t{geneend}\t{-used_strand}\t{revspecseq}\n")

        kmer_stroi.write(memchunk.getvalue())
         
        yield idx, cluster_dict, clusterpresab
 

def pattern_hasher(cluster_dict_iter, hash_pat, kmer_hash, genepres, patfilt, maf, output):
    #iterates through the cluster dictionary output by cluster_cutter()
    #outputs the hashed patterns, patterns and k-mers to the output files
    #two files are created: hashed k-mer patterns to presence/absence patterns
    #                       k-mer sequences to hashed k-mer patterns
    
    # set flag for multiple files
    multiple_files = False
    if hash_pat is None or kmer_hash is None:
        multiple_files = True

    if not multiple_files:
        memchonkheader1 = StringIO()
        memchonkheader2 = StringIO()
        
        memchonkheader1.write("hashed_pattern")
        for strain in sorted(genepres.columns):
            memchonkheader1.write(f"\t{strain}")
        memchonkheader1.write("\n")
        hash_pat.write(memchonkheader1.getvalue())

        memchonkheader2.write("k-mer\thashed_pattern\n")
        kmer_hash.write(memchonkheader2.getvalue())
    
    # keep track of already observed patterns
    # might have a big memory footprint
    # TODO: check for memory footprint
    patterns = set()

    for idx, cluster_dict, clusterpresab in cluster_dict_iter:
        if multiple_files:
            # close existing file?
            if hash_pat is not None:
                hash_pat.close()
            if kmer_hash is not None:
                kmer_hash.close()
            # we need to create the file
            path = os.path.join(output, idx)
            if not os.path.exists(path):
                os.mkdir(path)
            hash_pat, kmer_hash = create_hash_files(path)
            # must "flush" the already observed patterns file
            patterns = set()

            memchonkheader1 = StringIO()
            memchonkheader2 = StringIO()
            
            memchonkheader1.write("hashed_pattern")
            for strain in sorted(genepres.columns):
                memchonkheader1.write(f"\t{strain}")
            memchonkheader1.write("\n")
            hash_pat.write(memchonkheader1.getvalue())

            memchonkheader2.write("k-mer\thashed_pattern\n")
            kmer_hash.write(memchonkheader2.getvalue())

        memchunkhash_pat = StringIO()
        memchunkkmer_hash = StringIO()
        
        for kmer in cluster_dict:
            af = cluster_dict[kmer].sum() / cluster_dict[kmer].shape[0]
            if af >= 0.5:
                af = 1-af
            if af < maf:
                continue
            
            if patfilt == True:
                if tuple(cluster_dict[kmer]) == tuple(clusterpresab):
                    continue

            pattern = cluster_dict[kmer].view(np.uint8)
            hashed = hashlib.md5(pattern)
            khash = binascii.b2a_base64(hashed.digest()).decode()[:24]
            memchunkkmer_hash.write(f"{kmer}\t{khash}\n")
            
            
            if khash in patterns:
                continue
            patterns.add(khash)

            if not multiple_files and not len(patterns) % 1000:
                logger.debug(f"Observed patterns so far: {len(patterns)}")

            patterntup = "\t".join(map(str, cluster_dict[kmer]))
                        
            memchunkhash_pat.write(f"{khash}\t{patterntup}\n")

        hash_pat.write(memchunkhash_pat.getvalue())
        kmer_hash.write(memchunkkmer_hash.getvalue())


if __name__ == "__main__":
    pass
