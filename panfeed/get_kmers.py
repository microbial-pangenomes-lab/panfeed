#!/usr/bin/env python


import sys
import logging
import argparse
import numpy as np
import pandas as pd

from .__init__ import __version__
from .colorlog import ColorFormatter

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
    description = "Annotate association results with positional information"

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-a", "--associations",
            required=True,
            help="TSV file containing hashes and their significance "
                 "(e.g. pyseer output; tab-delimited, should have a header"
                 ", first column contains "
                 "the hash, and another column - by default 'lrt-pvalue' - "
                 "contains the association p-value, "
                 "which can be used to filter the table)")

    parser.add_argument("-p", "--kmers-to-hashes",
            required=True,
            help="TSV file relating gene clusters, kmers, and their hashes "
                 "(i.e. panfeed's kmers_to_hases.tsv file)")

    parser.add_argument("-k", "--kmers",
            required=True,
            help="TSV file with positional information of individual k-mers "
                 "(i.e. panfeed's kmers.tsv file)")

    parser.add_argument("-t", "--threshold",
            type=float,
            default=1,
            help="Association p-value threshold (default %(default).2f)")

    parser.add_argument("-c", "--column",
            default='lrt-pvalue',
            help="P-value column in the associations file (default %(default)s)")

    parser.add_argument("-o", "--output",
            default=None,
            help="Filename to save filtered associations table (not saved by default)")

    parser.add_argument("--only-passing",
            action="store_true",
            default=False,
            help="Only output passing k-mers (default is all)")

    parser.add_argument("-v", action='count',
                        default=0,
                        help='Increase verbosity level')
    parser.add_argument('--version', action='version',
                        version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    args = get_options()

    set_logging(args.v)

    a = pd.read_csv(args.associations, sep='\t', index_col=0)
    a.index.name = 'hashed_pattern'
    if args.column not in a.columns:
        logger.warning(f"Associations file does not have the {args.column} column")
        sys.exit(1)

    # filter associations table
    a = a[a[args.column] <= args.threshold]
    passing_hashes = set(a.index)

    logger.info(f"{len(passing_hashes)} patterns pass the p-value threshold {args.threshold}")
    if args.output is not None:
        a.to_csv(args.output, sep='\t')
        logger.info(f"Saved filtered associations to {args.output}")

    # load kmers to hashes table, piece-wise to reduce memory footprint
    # thanks to SO: https://stackoverflow.com/a/13653490/1237531
    iter_h = pd.read_csv(args.kmers_to_hashes, sep='\t',
                         iterator=True, chunksize=100_000)
    h = pd.concat([x[x['hashed_pattern'].isin(passing_hashes)]
        for x in iter_h]).set_index('hashed_pattern')
    clusters = set(h['cluster'].unique())
    kmers = set(h['k-mer'].unique())

    logger.info(f"Found {len(clusters)} gene clusters")
    logger.info(f"Found {len(kmers)} k-mers")

    # load k-mers table piece-wise to reduce memory footprint
    # thanks to SO: https://stackoverflow.com/a/13653490/1237531
    iter_k = pd.read_csv(args.kmers, sep='\t',
                         iterator=True, chunksize=100_000)
    k = pd.concat([x[x['cluster'].isin(clusters)]
        for x in iter_k]).set_index(['cluster', 'k-mer'])

    a = a.join(h, how='inner')
    if args.only_passing:
        how = 'left'
    else:
        how = 'right'
    a = a.reset_index().set_index(['cluster', 'k-mer']).join(k, how=how)

    a.to_csv(sys.stdout, sep='\t')


if __name__ == "__main__":
    main()
