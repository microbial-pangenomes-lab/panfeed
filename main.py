import pandas as pd
from collections import namedtuple
from pyfaidx import Fasta, Faidx
import numpy as np
import os, sys
import hashlib
import binascii
import itertools
from marcos_modules import create_faidx, parse_gff, iter_gene_clusters, Feature
from time import time
from new_modules import prep_data_n_fasta, what_are_my_inputfiles, set_input_output, cluster_cutter, pattern_hasher
import argparse

parser = argparse.ArgumentParser(description="input hybrid GFF/fasta files (prokka output). The output are 3 files that contain information for the specified strains of interest, hashed presence absence patterns to presence absence patterns and k-mers to presence absence patterns")

parser.add_argument("--stroi_in", type = str,
                    default="",
                    help = 'new line character delimited .txt file that contains the strain names, for which k-mer positions should be logged')

parser.add_argument("--stroi_out", type = str,
                    default = "k-mers_for_strains_of_interest.out",
                    help = "path and desired file name of output file to store k-mer information")

parser.add_argument("--kmer_length", type = int,
                    default = 31,
                    help = "k-mer length...")

parser.add_argument("--presence_absence", type = str,
                    default = "3gx1000gens.csv",
                    help = "the presence absence table as output by panaroo")

parser.add_argument("--gene_data", type = str,
                    default = "",
                    help = "gene data table as output by panaroo")

parser.add_argument("--start_inter", type = int,
                    default = 0,
                    help = "Interval to include upstream of the actual gene sequence (e.g. to include promoter region). default = 0")

parser.add_argument("--end_inter", type = int,
                    default = 0,
                    help = "Interval to include downstream of the actual gene sequence. default = 0")

args = parser.parse_args()


if __name__ == "__main__":

    strt = time()

    # Feature = namedtuple('Feature', ['id', 
    #                              'chromosome',
    #                              'start', 
    #                              'end', 
    #                              'strand'])

    klength = args.kmer_length

    filelist = what_are_my_inputfiles()

    data = {}

    stroi, error_log, kmer_stroi, hash_pat, kmer_hash, genedata, genepres = set_input_output(args.stroi_in, 
                                                                                             args.stroi_out, 
                                                                                             args.presence_absence, 
                                                                                             args.gene_data)

    sys.stderr = error_log

    data = prep_data_n_fasta(filelist)
    
    kmer_pattern_dict = {}
    


    pattern_hasher(cluster_cutter(iter_gene_clusters(genepres, 
                                                     data,
                                                     error_log,
                                                     args.start_inter,
                                                     args.end_inter),
                                  klength, 
                                  stroi, 
                                  kmer_stroi), 
                   hash_pat, 
                   kmer_hash,
                   genepres)




# Parallel(n_jobs=2, backend="threading", batch_size = 1)(delayed(pattern_hasher)(cluster_dict, hash_pat, kmer_hash) for cluster_dict in Parallel(n_jobs=2, backend="threading", batch_size = 1)(delayed(cluster_cutter)(gene_cluster, klength, stroi, kmer_stroi) for gene_cluster in iter_gene_clusters(genepres, data)))

    nd = time()

    print(nd - strt)

    kmer_stroi.close()

    hash_pat.close()

    kmer_hash.close()

    # stroi.close()


#3gxallgen = 241s   229s
#2gxallgen = 155s   165s
#1gxallgen = 98s    94s































