#!/usr/bin/env python


import sys
import logging
import argparse
import logging.handlers

from .__init__ import __version__
from .colorlog import ColorFormatter

from .panfeed import prep_data_n_fasta, what_are_my_inputfiles, set_input_output
from .panfeed import clean_up_fasta, what_are_my_inputfiles, set_input_output
from .panfeed import cluster_cutter, pattern_hasher, create_faidx, parse_gff
from .panfeed import iter_gene_clusters, Feature


logger = logging.getLogger('panfeed')


def set_logging(v):
    logger.propagate = True
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    if v == 0:
        ch.setLevel(logging.INFO)
    elif v >= 1:
        ch.setLevel(logging.DEBUG)
    formatter = ColorFormatter('%(asctime)s - %(name)s - $COLOR%(message)s$RESET','%H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


def get_options():
    description = "Get gene cluster specific k-mers from a set of bacterial genomes"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-g", "--gff",
                        help = "Directory containing all samples' GFF "
                               "files (must contain nucleotide sequence as "
                               "well, and samples should be named in the "
                               "same way as in the panaroo header)")

    parser.add_argument("--targets",
                        default=None,
                        help = "File indicating for which samples "
                               "k-mer positions should be logged "
                               "(one sample name per line, "
                               "default: no samples)")

    parser.add_argument("-o", "--output",
                        default = "panfeed",
                        help = "Output directory to store outputs "
                               "(will cause an error if already present)")

    parser.add_argument("-k", "--kmer-length", type = int,
                        default = 31,
                        help = "K-mer length (default: %(default)d)")

    parser.add_argument("-p", "--presence-absence",
                        help = "Gene clusters presence absence table "
                               "as output by panaroo")

    parser.add_argument("--maf", type = float,
                        default = 0.01,
                        help = "Minor allele frequency threshold; "
                               "patterns whose frequency is below "
                               "this value or above 1-MAF are excluded "
                               "(default: %(default).2f, does not apply "
                               "to the kmers.tsv file)")
    
    parser.add_argument("--upstream", type = int,
                        default = 0,
                        help = "How many bases to include upstream of "
                               "the actual gene sequences "
                               "(e.g. to include the 5' region, "
                               "default: %(default)d)")
    
    parser.add_argument("--downstream", type = int,
                        default = 0,
                        help = "How many bases to include downstream of "
                               "the actual gene sequences "
                               "(e.g. to include the 3' region, "
                               "default: %(default)d)")

    parser.add_argument("--downstream-start-codon",
                        action = "store_true",
                        default = False,
                        help = "Center the --downstream argument "
                               "at the start codon (default is stop codon)")

    parser.add_argument("--non-canonical",
                        action = "store_true",
                        default = False,
                        help = "Compute non-canonical k-mers "
                               "(default is canonical)")

    parser.add_argument("--no-filter",
                        action = "store_true",
                        default = False,
                        help = "Do not filter out k-mers with the same "
                               "presence absence pattern as the gene "
                               "cluster itself")

    parser.add_argument("--multiple-files",
                        action = "store_true",
                        default = False,
                        help = "Generate one set of outputs for each "
                               "gene cluster (default: one set of outputs)")

    parser.add_argument("-v", action='count',
                        default=0,
                        help='Increase verbosity level')
    parser.add_argument('--version', action='version',
                        version='%(prog)s '+__version__)
    
    return parser.parse_args()


def main():
    args = get_options()

    set_logging(args.v)    

    klength = args.kmer_length
    
    if args.downstream_start_codon == True and args.upstream + args.downstream < klength:
        logger.warning("Query sequence is shorter than k-mer length"
                       "Decrease k-mer size or increase query sequence length")
        sys.exit(1)
    
    if args.maf > 0.5:
        logger.warning("--maf should be below 0.5")
        sys.exit(1)

    logger.info("Looking at input GFF files")
    filelist = what_are_my_inputfiles(args.gff)
    logger.info(f"Found {len(filelist)} input genomes")

    data = {}

    logger.info("Preparing output files")
    (stroi, kmer_stroi, hash_pat,
     kmer_hash, genepres) = set_input_output(args.targets,
                                             args.presence_absence,
                                             args.output,
                                             not args.multiple_files)

    logger.info("Preparing inputs")
    data = prep_data_n_fasta(filelist, args.gff, args.output)
    
    kmer_pattern_dict = {}
    
    logger.info("Extracting k-mers")
    pattern_hasher(cluster_cutter(iter_gene_clusters(genepres, 
                                                     data,
                                                     args.upstream,
                                                     args.downstream,
                                                     args.downstream_start_codon,
                                                     not args.no_filter),
                                  klength,
                                  stroi, 
                                  kmer_stroi,
                                  not args.non_canonical,
                                  args.output),
                   hash_pat, 
                   kmer_hash,
                   genepres,
                   not args.no_filter,
                   args.maf,
                   args.output)

    # nd = time()

    # print(nd - strt)

    if kmer_stroi is not None:
        kmer_stroi.close()

    if hash_pat is not None:
        hash_pat.close()

    if kmer_hash is not None:
        kmer_hash.close()
    
    logger.info("Removing temporary fasta files and faidx indices")
    clean_up_fasta(filelist, args.output)


if __name__ == "__main__":
    main()
