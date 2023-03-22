#!/usr/bin/env python


import sys
import logging
import argparse
import logging.handlers
from functools import partial
import itertools
from multiprocessing import Process, Queue

from .__init__ import __version__
from .colorlog import ColorFormatter

from .input import prep_data_n_fasta, what_are_my_inputfiles, set_input_output
from .input import clean_up_fasta, what_are_my_inputfiles, set_input_output
from .panfeed import cluster_cutter, pattern_hasher, write_headers
from .input import iter_gene_clusters


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


def worker(f, read_q, write_q):
    while True:
        work = read_q.get()
        if work is None:
            logger.debug("Worker process finished")
            write_q.put(None)
            return

        ret = f(work)
        if len(ret) == 0:
            logger.debug("Worker process finished")
            write_q.put(None)
            return
        write_q.put((ret,),timeout = 300)


def reader(iter_i, read_q, n):
    for x in iter_i:
        read_q.put(x)

    logger.debug("Reader process finished")

    # poison pill
    # one for each worker
    for _ in range(n):
        read_q.put(None)


def writer(f, write_q, n):
    dead = 0

    patterns = set()

    while True:
        work = write_q.get()
        if work is None:
            dead += 1
            logger.debug(f"Finished worker processes {dead}/{n}")
            if dead == n:
                logger.debug("All worker processes are finished")
                return
            continue
        patterns = f(work, patterns=patterns)


def get_options():
    description = "Get gene cluster specific k-mers from a set of bacterial genomes"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-g", "--gff",
                        required=True,
                        help = "Directory containing all samples' GFF "
                               "files (must contain nucleotide sequence as "
                               "well unless -f is used, "
                               "and samples should be named in the "
                               "same way as in the panaroo header)")

    parser.add_argument("-p", "--presence-absence",
                        required=True,
                        help = "Gene clusters presence absence table "
                               "as output by panaroo")

    parser.add_argument("--targets",
                        default=None,
                        help = "File indicating for which samples "
                               "k-mer positions should be logged "
                               "(one sample name per line, "
                               "default: no samples)")

    parser.add_argument("--genes",
                        default=None,
                        help = "File indicating for which gene clusters "
                               "k-mer positions and patterns should be logged "
                               "(one gene cluster ID per line, "
                               "default: all gene clusters)")

    parser.add_argument("-o", "--output",
                        default = "panfeed",
                        help = "Output directory to store outputs "
                               "(will cause an error if already present)")

    parser.add_argument("-f", "--fasta",
                        help = "Directory containing all samples' nucleotide "
                               "fasta files (extension either .fasta "
                               "or .fna, "
                               "samples should be named in the "
                               "same way as in the panaroo header")

    parser.add_argument("-k", "--kmer-length", type = int,
                        default = 31,
                        help = "K-mer length (default: %(default)d)")

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

    parser.add_argument("--consider-missing",
                        action = "store_true",
                        default = False,
                        help = "Output NaN for strains that do not encode "
                               "for a k-mer if the gene is missing "
                               "(default: value is set to 0, as for the gene. "
                               "WARNING: slows down execution)")

    parser.add_argument("--multiple-files",
                        action = "store_true",
                        default = False,
                        help = "Generate one set of outputs for each "
                               "gene cluster (default: one set of outputs)")

    parser.add_argument("--compress",
                        action = "store_true",
                        default = False,
                        help = "Compress output files with gzip "
                               "(default: plain text)")

    parser.add_argument("--cores",
                        type = int,
                        default = 1,
                        help = "Threads (default: %(default)d, at least 3 "
                               "are needed for parallelization)")

    parser.add_argument("-ql", "--queue-limit",
                        type = int,
                        default = 3,
                        help = "limit on items that may be put into the reading"
                               "and writing queues. (default: %(default)d)"
                               "this option is only relevant for cores > 1"
                               "reading queue limit = ql * cores"
                               "writing queue limit = ql")

    parser.add_argument("-v", action='count',
                        default=0,
                        help='Increase verbosity level')
    parser.add_argument('--version', action='version',
                        version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    args = get_options()

    set_logging(args.v)

    qlimit = args.queue_limit

    klength = args.kmer_length

    if args.downstream_start_codon == True and args.upstream + args.downstream < klength:
        logger.warning("Query sequence is shorter than k-mer length"
                       "Decrease k-mer size or increase query sequence length")
        sys.exit(1)

    if args.maf > 0.5:
        logger.warning("--maf should be below 0.5")
        sys.exit(1)

    logger.info("Looking at input GFF files")
    filelist, fastalist = what_are_my_inputfiles(args.gff, args.fasta)
    logger.info(f"Found {len(filelist)} input genomes")
    if args.fasta is not None:
        logger.info(f"Found {len(fastalist)} input nucleotide sequences")

    data = {}

    logger.info("Preparing output files")
    (stroi, genes, kmer_stroi, hash_pat,
     kmer_hash, genepres) = set_input_output(args.targets,
                                             args.genes,
                                             args.presence_absence,
                                             args.output,
                                             not args.multiple_files,
                                             args.compress)

    logger.info("Preparing inputs")
    data = prep_data_n_fasta(filelist, fastalist,
                             args.gff, args.fasta, args.output)

    if not args.multiple_files:
        write_headers(hash_pat, kmer_hash, genepres)

    logger.info("Extracting k-mers")
    iter_i = iter_gene_clusters(genepres, 
                                data,
                                args.upstream,
                                args.downstream,
                                args.downstream_start_codon,
                                not args.no_filter,
                                genes)
    iter_o = partial(cluster_cutter,
                     klength=klength,
                     stroi=stroi,
                     multiple_files=args.multiple_files,
                     canon=not args.non_canonical,
                     consider_missing_cluster=args.consider_missing,
                     output=args.output,
                     compress=args.compress)

    patterns = set()
    func_w = partial(pattern_hasher,
                     kmer_stroi=kmer_stroi,
                     hash_pat=hash_pat,
                     kmer_hash=kmer_hash,
                     genepres=genepres,
                     patfilt=not args.no_filter,
                     maf=args.maf,
                     consider_missing_cluster=args.consider_missing,
                     output=args.output,
                     patterns=patterns,
                     compress=args.compress)

    if args.cores > 2:
        # thanks to @SamStudio8 for the inspiration
        read_q = Queue(maxsize = qlimit * args.cores)
        write_q = Queue(maxsize = qlimit)

        processes = []

        reader_process = Process(
            target=reader,
            args=(
                iter_i,
                read_q,
                args.cores - 2
            ),
        )
        processes.append(reader_process)

        writer_process = Process(
            target=writer,
            args=(
                func_w,
                write_q,
                args.cores - 2,
            ),
        )
        processes.append(writer_process)

        for i in range(args.cores - 2):
            p = Process(
                target=worker,
                args=(
                    iter_o,
                    read_q,
                    write_q,
                ),
            )
            processes.append(p)
        logger.debug(f"Started {args.cores - 2} worker processes")

        for p in processes:
            p.start()

        for p in processes:
            p.join()
    else:
        if args.cores > 1:
            logger.warning("Need at least 3 cores for parallelization")
            logger.warning("Running a single thread")

        patterns = set()
        for x in iter_i:
            ret = iter_o(x)
            if len(ret) == 0:
                continue
            patterns = func_w((ret,),
                              patterns=patterns)


    if kmer_stroi is not None:
        kmer_stroi.close()

    if hash_pat is not None:
        hash_pat.close()

    if kmer_hash is not None:
        kmer_hash.close()

    logger.info("Removing temporary fasta files and faidx indices")
    clean_up_fasta(filelist, fastalist, args.output, args.fasta)


if __name__ == "__main__":
    main()
