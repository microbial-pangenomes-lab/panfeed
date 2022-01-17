import pandas as pd
from collections import namedtuple
from pyfaidx import Fasta, Faidx
import numpy as np
import os, sys
import hashlib
import binascii
import itertools
# from time import time
from panfeed_modules import prep_data_n_fasta, what_are_my_inputfiles, set_input_output, cluster_cutter, pattern_hasher, create_faidx, parse_gff, iter_gene_clusters, Feature
import argparse

parser = argparse.ArgumentParser(description="input hybrid GFF/fasta files (prokka output). The output are 3 files that contain information for the specified strains of interest, hashed presence absence patterns to presence absence patterns and k-mers to presence absence patterns")

parser.add_argument("--stroi_in", type = str,
                    default="stroi.txt",
                    help = 'new line character delimited .txt file that contains the strain names, for which k-mer positions should be logged')

parser.add_argument("--stroi_out", type = str,
                    default = "k-mers_for_strains_of_interest.csv",
                    help = "path and desired file name of output file to store k-mer information")

parser.add_argument("--kmer_length", type = int,
                    default = 31,
                    help = "k-mer length...")

parser.add_argument("--presence_absence", type = str,
                    default = "3gx10gens.csv",
                    help = "the presence absence table as output by panaroo")

parser.add_argument("--gene_data", type = str,
                    default = "",
                    help = "gene data table as output by panaroo")

parser.add_argument("--start_inter", type = int,
                    default = 0,
                    help = "Interval to include upstream of the actual gene sequence (e.g. to include promoter region). default = 0")

parser.add_argument("--end_inter", type = int,
                    default = 0,
                    help = "Interval to include downstream of the actual gene sequence. Default = 0")

parser.add_argument("--canon", type = bool,
                    default = True,
                    help = "Switching between counting canonical (True) and non-canonical k-mers (False). Default is True")

parser.add_argument("--specfilt", type = bool,
                    default = True,
                    help = "If True: filters out all kmers, from one gene cluster, with the same presence absence pattern as the gene cluster itself. Default is True")

args = parser.parse_args()


if __name__ == "__main__":

    # strt = time()

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
                                                     args.end_inter,
                                                     args.canon,
                                                     args.specfilt),
                                  klength, 
                                  stroi, 
                                  kmer_stroi,
                                  args.canon), 
                   hash_pat, 
                   kmer_hash,
                   genepres,
                   args.specfilt)


    # nd = time()

    # print(nd - strt)

    kmer_stroi.close()

    hash_pat.close()

    kmer_hash.close()





























