#!/usr/bin/env python


import sys
import logging
import argparse
import numpy as np
import pandas as pd

from matplotlib import colors
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
    description = "Plot association results from panfeed"

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-k", "--kmers",
            required=True,
            help="TSV file containing the output of panfeed-get-kmers")

    parser.add_argument("-c", "--column",
            default="lrt-pvalue",
            help="P-value column in the associations file (default %(default)s)")

    parser.add_argument("-t", "--threshold",
            type=float,
            default=1,
            help="Association p-value threshold (default %(default).2f)")

    parser.add_argument("-p", "--phenotype",
            required=True,
            help="Phenotype file in TSV format, used to list all strains")

    parser.add_argument("--phenotype-column",
            default=None,
            help="Column in phenotype TSV file for sorting "
                 "(default is sorting by p-values)")

    parser.add_argument('--sample',
                        type=float,
                        default=None,
                        help="Only show a randomly picked set of sample "
                             "(a value between 0 and 1 indicating the proportion"
                             " to show, default all)")

    parser.add_argument('--start',
                        type=int,
                        default=None,
                        help="Relative position to start the plots "
                             "(default all available positions)")

    parser.add_argument('--stop',
                        type=int,
                        default=None,
                        help="Relative position to end the plots "
                             "(default all available positions)")

    parser.add_argument('--format',
                        choices=('png',
                                 'tiff',
                                 'pdf',
                                 'svg'),
                        default='png',
                        help='Output format for plots (default %(default)s)')

    parser.add_argument('--dpi',
                        type=int,
                        default=300,
                        help='Output resolution (DPI, default %(default)d)')

    parser.add_argument("--minimum-pvalue",
            type=float,
            default=1E-10,
            help="Minimum p-value for color and transparency (default %(default).2e)")

    parser.add_argument("--nucleotides",
            action="store_true",
            default=False,
            help="Draw nucleotide sequence on all plots "
                 "(WARNING: only makes sense if few samples and positions are considered)")

    parser.add_argument('--alpha',
                        type=float,
                        default=0,
                        help="Opacity for non-passing k-mers "
                             "(between 0 and 1, 0 indicates full transparency"
                             ", default %(default).2f)")

    parser.add_argument('--xticks',
                        type=int,
                        default=200,
                        help="Spacing for ticks on x axis (default %(default)d)")

    parser.add_argument('--height',
                        type=int,
                        default=9,
                        help="Figure height (inches, default %(default)d)")

    parser.add_argument('--width',
                        type=int,
                        default=10,
                        help="Figure width (inches, default %(default)d)")

    parser.add_argument("-v", action='count',
                        default=0,
                        help='Increase verbosity level')
    parser.add_argument('--version', action='version',
                        version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    args = get_options()

    set_logging(args.v)

    if args.sample is not None and (args.sample > 1 or args.sample < 0):
        logger.warning("--sample should be between 0 and 1")
        sys.exit(1)
    if args.alpha > 1 or args.alpha < 0:
        logger.warning("--alpha should be between 0 and 1")
        sys.exit(1)
    if (args.start is not  None and args.stop is None) or (args.start is None and args.stop is not None):
        logger.warning("both --start and --stop are needed")
        sys.exit(1)
    if args.start is not None and args.start > args.stop:
        logger.warning("--start should be lower than --stop")
        sys.exit(1)
    if args.nucleotides and (args.sample is None or args.start is None):
        logger.warning("drawing nucleotide sequences without zooming in might "
                       "increase plotting time and memory consumption "
                       "while generating useless plots")

    p = pd.read_csv(args.phenotype, sep='\t', index_col=0)
    pbinary = False
    if args.phenotype_column is not None:
        if args.phenotype_column not in p.columns:
            logger.warning(f"phenotype file does not have the {args.phenotype_column} column")
            sys.exit(1)
        p = p[args.phenotype_column].dropna().sort_values(ascending=False)
        if args.sample is not None:
            p = p.sample(frac=args.sample)
        # guess is binary or continuous
        pvalues = set(p.values)
        if len(pvalues) == 2 and 1 in pvalues and 0  in pvalues:
            logger.info("Phenotype is binary")
            pbinary = True
            # index for transition from 1 to 0
            yindex = [i for i, x in enumerate(p)
                      if i > 0
                      and x != p.iloc[i - 1]]
            if len(yindex) > 0:
                yindex = yindex[0]
            else:
                yindex = None
        else:
            logger.info("Phenotype is continuos")
        logger.info(f"sorting samples by their {args.phenotype_column} phenotype")
    else:
        logger.info(f"sorting samples by their lowest association p-value")
        if args.sample is not None:
            p = p.sample(frac=args.sample)
    strains = set(p.index)
    logger.info(f"plots will include {len(strains)} samples")

    k = pd.read_csv(args.kmers, sep='\t')
    if args.column not in k.columns:
        logger.warning(f"k-mer file does not have the {args.column} column")
        sys.exit(1)
    # ignore samples not in phenotype file
    k = k[k['strain'].isin(strains)]
    # compute -log10 of p-value
    k['significance'] = -np.log10(k[args.column])
    clusters = set(k['cluster'].unique())

    bases = {'A': 'T',
             'T': 'A',
             'G': 'C',
             'C': 'G'}
    base2int = {'A': 0,
                'G': 1,
                'T': 2,
                'C': 3}

    # add a column with the nucleotide
    k['base'] = np.nan
    k['base'] = k['k-mer'].str.upper().str[0]
    k.loc[k[k['strand'] == -1].index,
          'base'] = [bases.get(x[-1], 'N')
                     for x in k[k['strand'] == -1
                         ]['k-mer'].str.upper().values]
    # convert it to scalar
    k['scalar'] = [base2int.get(x, np.nan) for x in k['base'].values]

    logger.info(f"preparing plots for {len(clusters)} gene clusters")

    def handle_paralogs(x):
        if len(x) > 1:
            return 99
        return x

    def handle_paralogs_text(x):
        if len(x) > 1:
            return '-'
        return x

    cmap1 = plt.get_cmap('viridis').copy()
    cmap1.set_bad('xkcd:grey')
    cmap1.set_under('xkcd:light grey')

    base_colors = list(sns.color_palette('tab20', 4))
    cmap2 = colors.LinearSegmentedColormap.from_list('nucleotides', base_colors, 4)
    cmap2.set_bad('xkcd:grey')
    cmap2.set_over(list(sns.color_palette('tab20', 5))[-1])

    for gene in sorted(clusters):
        logger.info(f"preparing plots for {gene}")

        if args.start is not None:
            # zoom in
            cl = k[(k['cluster'] == gene) &
                   (k['gene_start'] >= args.start) &
                   (k['gene_start'] <= args.stop)]
        else:
            cl = k[k['cluster'] == gene]
        # reality check
        if cl.shape[0] == 0:
            logger.warning(f"Skipping {gene}")
            continue

        # p-value matrix
        g = cl.pivot_table(index='strain', columns='gene_start',
                           values='significance',
                           aggfunc=max)
        # add missing strains for which no entry is present in k-mers table
        if len(strains.difference(g.index)) > 0:
            g = g.reindex(sorted(g.index) + sorted(strains.difference(g.index)))
        if args.phenotype_column is None:
            # sort the table by p-values
            g = g.loc[g.fillna(0).T.max().sort_values(ascending=False).index]
        else:
            # sort by phenotype
            g = g.loc[p.index]

        # nucleotide matrix
        b = cl.pivot_table(index='strain', columns='gene_start',
                           values='scalar',
                           aggfunc=handle_paralogs)
        # add missing strains for which no entry is present in k-mers table
        if len(strains.difference(b.index)) > 0:
            b = b.reindex(sorted(b.index) + sorted(strains.difference(b.index)))
        if args.phenotype_column is None:
            # sort the table by p-values
            b = b.loc[g.index]
        else:
            # sort by phenotype
            b = b.loc[p.index]
        if args.nucleotides:
            # nucleotide matrix
            t = cl.pivot_table(index='strain', columns='gene_start',
                               values='base',
                               aggfunc=handle_paralogs_text)
            # add missing strains for which no entry is present in k-mers table
            if len(strains.difference(t.index)) > 0:
                t = t.reindex(sorted(t.index) + sorted(strains.difference(t.index)),
                              fill_value='-')
            if args.phenotype_column is None:
                # sort the table by p-values
                t = t.loc[g.index]
            else:
                # sort by phenotype
                t = t.loc[p.index]

        # find the index for the gene start (if present)
        gene_start = [i for i, x in enumerate(g.columns)
                      if x == 0]
        if len(gene_start) > 0:
            gene_start = gene_start[0]
        else:
            gene_start = None
        # create an index to draw xticks as relative gene positions
        xticks = [(i, x) for i, x in enumerate(g.columns)
                  if not x % args.xticks]

        fig, ax = plt.subplots(figsize=(args.width, args.height))
        im = ax.imshow(g.values, cmap=cmap1, vmin=-np.log10(args.threshold),
                       vmax=-np.log10(args.minimum_pvalue),
                       aspect="auto", interpolation='none',
                       rasterized=True, alpha=1)
        if args.nucleotides:
            for x in range(t.shape[1]):
                for y in range(t.shape[0]):
                    ax.text(x, y, t.iloc[y, x],
                            ha='center', va='center')
        ax.set_yticks([])
        if pbinary and yindex is not None:
            ax.axhline(yindex, lw=1, color="black")
        if gene_start is not None:
            ax.axvline(gene_start, lw=1, color="black")
        ax.set_xticks([x[0] for x in xticks],
                      labels=[x[1] for x in xticks])
        ax.set_title(f"significant k-mers {gene}")
        ax.set_xlabel("position relative to gene start")
        ax.set_ylabel(f"{g.shape[0]} samples")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2.5%", pad=0.05)
        coba = plt.colorbar(im, cax=cax)
        coba.set_label("-log10 p-value")
        plt.savefig(f'significance_{gene}.{args.format}', dpi=args.dpi, bbox_inches="tight")
        logger.info(f'saved plot significance_{gene}.{args.format}')
        plt.close()

        fig, ax = plt.subplots(figsize=(args.width, args.height))
        im = ax.imshow(b.values, cmap=cmap2, vmin=0, vmax=3,
                       aspect="auto", interpolation='none',
                       rasterized=True, alpha=1)
        if args.nucleotides:
            for x in range(t.shape[1]):
                for y in range(t.shape[0]):
                    ax.text(x, y, t.iloc[y, x],
                            ha='center', va='center')
        ax.set_yticks([])
        if pbinary and yindex is not None:
            ax.axhline(yindex, lw=1, color="black")
        if gene_start is not None:
            ax.axvline(gene_start, lw=1, color="black")
        ax.set_xticks([x[0] for x in xticks],
                      labels=[x[1] for x in xticks])
        ax.set_title(f"nucleotide sequence {gene}")
        ax.set_xlabel("position relative to gene start")
        ax.set_ylabel(f"{g.shape[0]} samples")
        plt.savefig(f'sequence_{gene}.{args.format}', dpi=args.dpi, bbox_inches="tight")
        logger.info(f'saved plot sequence_{gene}.{args.format}')
        plt.close()

        g = g.fillna(g.min().min())
        g = (g - -np.log10(args.threshold)) / (-np.log10(args.minimum_pvalue) - -np.log10(args.threshold))
        g[g < args.alpha] = args.alpha
        g[g > 1] = 1
        g[np.isnan(g)] = args.alpha

        fig, ax = plt.subplots(figsize=(args.width, args.height))
        im = ax.imshow(b.values, cmap=cmap2, vmin=0, vmax=3,
                       aspect="auto", interpolation='none',
                       rasterized=True, alpha=g.values)
        if args.nucleotides:
            for x in range(t.shape[1]):
                for y in range(t.shape[0]):
                    ax.text(x, y, t.iloc[y, x],
                            ha='center', va='center')
        ax.set_yticks([])
        if pbinary and yindex is not None:
            ax.axhline(yindex, lw=1, color="black")
        if gene_start is not None:
            ax.axvline(gene_start, lw=1, color="black")
        ax.set_xticks([x[0] for x in xticks],
                      labels=[x[1] for x in xticks])
        ax.set_title(f"significant k-mers {gene}")
        ax.set_xlabel("position relative to gene start")
        ax.set_ylabel(f"{g.shape[0]} samples")
        plt.savefig(f'hybrid_{gene}.{args.format}', dpi=args.dpi, bbox_inches="tight")
        logger.info(f'saved plot hybrid_{gene}.{args.format}')
        plt.close()

    # plot nucleotide color scheme
    plt.figure(figsize=(3, 0.7))

    sns.heatmap(data=[[0, 1, 2, 3]],
                annot=[['A', 'G', 'T', 'C']],
                fmt='s',
                annot_kws={'size': 16, 'weight': 'bold'},
                cmap=cmap2,
                cbar=False,
                yticklabels=False,
                xticklabels=False,
                square=True)

    plt.savefig(f'sequence_legend.{args.format}', dpi=args.dpi, bbox_inches="tight")
    logger.info(f"Saved nucleotide color key at sequence_legend.{args.format}")
    plt.close()


if __name__ == "__main__":
    main()
