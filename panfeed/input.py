#!/usr/bin/env python

import os
import sys
import gzip
import logging
import numpy as np
import pandas as pd
from pyfaidx import Fasta

from .classes import Feature, Seqinfo

logger = logging.getLogger('panfeed.input')


def what_are_my_inputfiles(gffdir, fastadir=None): # adds all the GFF files in the directory to a list
    filelist = set()
    fastalist = set()

    if os.path.isfile(gffdir):
        # we have a file of files
        logger.debug(f"GFF input ({gffdir}) is a file of files")
        files = {x.rstrip() for x in open(gffdir)}
    else:
        # we have a directory
        logger.debug(f"GFF input ({gffdir}) is a directory")
        files = os.listdir(gffdir)
    for file in files:
        only_file = os.path.split(file)[-1]
        if file.endswith(".gff"):
            genome = ".".join(only_file.split(".")[:-1])
            logger.debug(f"Adding {genome} ({file})")
            filelist.add(genome)
        else:
            logger.debug(f"Ignoring {file}")

    if fastadir is not None:
        if os.path.isfile(fastadir):
            # we have a file of files
            logger.debug(f"Fasta input ({fastadir}) is a file of files")
            files = {x.rstrip() for x in open(fastadir)}
        else:
            # we have a directory
            logger.debug(f"Fasta input ({fastadir}) is a directory")
            files = os.listdir(fastadir)
        # check which of the GFF files have a corresponding
        # nucleotide fasta file
        for file in files:
            only_file = os.path.split(file)[-1]
            if file.endswith(".fasta") or file.endswith(".fna"):
                genome = ".".join(only_file.split(".")[:-1])
                if genome in filelist:
                    logger.debug(f"Adding {genome} nucleotides ({file})")
                    fastalist.add(genome)
                else:
                    logger.debug(f"Ignoring {genome} (no corresponding gff)")
            else:
                logger.debug(f"Ignoring {file}")


    if len(filelist) == 0:
        logger.error(f"No GFF files found in the inputs provided ({gffdir})")
        sys.exit(1)

    return sorted(filelist), sorted(fastalist)


def prep_data_n_fasta(filelist, fastalist,
                      gffdir, fastadir, output): # prepares the hybrid GFF files and creates fasta files
    data = {}

    if os.path.isfile(gffdir):
        # file of files
        # create a dict with the mapping between samples and full path
        gff_d = {}
        for file in open(gffdir):
            file = file.rstrip()
            only_file = os.path.split(file)[-1]
            if file.endswith(".gff"):
                genome = ".".join(only_file.split(".")[:-1])
                gff_d[genome] = file
    if fastadir is not None and os.path.isfile(fastadir):
        # file of files
        # create a dict with the mapping between samples and full path
        fasta_d = {}
        for file in open(fastadir):
            file = file.rstrip()
            only_file = os.path.split(file)[-1]
            if file.endswith(".fna") or file.endswith(".fasta"):
                genome = ".".join(only_file.split(".")[:-1])
                fasta_d[genome] = file

    for genome in filelist:
        logger.debug(f"Handling {genome}")

        if genome not in fastalist:
            logger.debug(f"Creating fasta file for {genome}")

            try:
                os.mkdir(os.path.join(output, 'fastas'))
            except FileExistsError:
                pass
            txtfile = open(os.path.join(output, 'fastas', f"{genome}.fasta"), "w")
            if os.path.isfile(gffdir):
                txtfile.write(open(gff_d[genome], "r").read().split("##FASTA")[1])
            else:
                txtfile.write(open(os.path.join(gffdir, f"{genome}.gff"),
                    "r").read().split("##FASTA")[1])
            txtfile.close()

            logger.debug(f"Creating faidx file for {genome}")
            sequences = create_faidx(os.path.join(output, 'fastas', f'{genome}.fasta'),
                                     close=True)
        else:
            if fastadir is not None and os.path.isfile(fastadir):
                ffile = fasta_d[genome]
            else:
                logger.debug(f"Creating faidx file for existing nucleotide sequence {genome}")
                file1 = os.path.join(fastadir, f'{genome}.fna')
                file2 = os.path.join(fastadir, f'{genome}.fasta')
                if os.path.exists(file1):
                    ffile = file1
                elif os.path.exists(file2):
                    ffile = file2
                else:
                    logger.error(f"Neither {genome}.fna not {genome}.fasta found in {fastadir}")
                    sys.exit(1)
            sequences = create_faidx(ffile, close=True)

        logger.debug(f"Parsing features for {genome}")
        if os.path.isfile(gffdir):
            features = parse_gff(gff_d[genome])
        else:
            features = parse_gff(os.path.join(gffdir, f'{genome}.gff'))

        data[genome] = (sequences, features)

    return data


def clean_up_fasta(filelist, fastalist, output, fastadir):
    if fastadir is not None and os.path.isfile(fastadir):
        # file of files
        # create a dict with the mapping between samples and full path
        fasta_d = {}
        for file in open(fastadir):
            file = file.rstrip()
            only_file = os.path.split(file)[-1]
            if file.endswith(".fna") or file.endswith(".fasta"):
                genome = ".".join(only_file.split(".")[:-1])
                fasta_d[genome] = file

    for genome in filelist:
        if genome in fastalist:
            logger.debug(f"Keeping fasta file for {genome}")
        else:
            fasta_file = os.path.join(output, 'fastas', f"{genome}.fasta")
            logger.debug(f"Removing fasta file for {genome} ({fasta_file})")
            if os.path.exists(fasta_file):
                os.remove(fasta_file)
            else:
                logger.warning(f"Could not delete {fasta_file}")

        if genome not in fastalist:
            faidx_file = os.path.join(output, 'fastas', f"{genome}.fasta.fai")
        else:
            if fastadir is not None and os.path.isfile(fastadir):
                faidx_file = fasta_d[genome] + ".fai"
            else:
                faidx_file1 = os.path.join(fastadir, 'fastas', f"{genome}.fasta.fai")
                faidx_file2 = os.path.join(fastadir, 'fastas', f"{genome}.fna.fai")
                if os.path.isfile(faidx_file1):
                    faidx_file = faidx_file1
                else:
                    faidx_file = faidx_file2
        logger.debug(f"Removing faidx file for {genome} ({faidx_file})")
        if os.path.exists(faidx_file):
            os.remove(faidx_file)
        else:
            logger.warning(f"Could not delete {faidx_file}")


def set_input_output(stroi_in, genes_in, presence_absence, output,
                     single_file=True, compress=False):
    # determines the names of the input, output files 

    logger.debug(f"Loading pangenome file from panaroo ({presence_absence})")
    genepres = pd.read_csv(presence_absence,
                           sep=",", index_col=0,
                           low_memory=False).drop(
                                   columns=['Non-unique Gene name', 'Annotation'])
    # load the gene_presence_absence matrix which will be used to ascribe strain/genes to clusters

    if stroi_in is not None:
        logger.debug(f"Loading target strains ({stroi_in})")
        stroi = set()
        strainfile = open(f"{stroi_in}", "r")
        for line in strainfile:
            stroi.add(line.rstrip("\n"))
    else:
        logger.warning(f"No target strains provided")
        stroi = ""

    if genes_in is not None:
        logger.debug(f"Loading target gene clusters ({genes_in})")
        genes = set()
        genesfile = open(f"{genes_in}", "r")
        for line in genesfile:
            genes.add(line.rstrip("\n"))
    else:
        genes = None

    if os.path.exists(output):
        logger.error(f"Output directory {output} exists! Please "
                      "remove it and restart")
        sys.exit(1)
    else:
        logger.debug(f"Creating output directory ({output})")
        os.mkdir(output)

    logger.debug(f"Creating output files within {output}")

    if single_file:
        kmer_stroi = create_kmer_stroi(output, compress)
        hash_pat, kmer_hash = create_hash_files(output, compress) 
    else:
        # delayied opening of the files
        kmer_stroi = None
        hash_pat = None
        kmer_hash = None

    return stroi, genes, kmer_stroi, hash_pat, kmer_hash, genepres


def create_kmer_stroi(output, compress=False):
    if not compress:
        kmer_stroi = open(os.path.join(output, "kmers.tsv"), "w")
    else:
        kmer_stroi = gzip.open(os.path.join(output, "kmers.tsv.gz"), "wt",
                               compresslevel=9)

    #creates the header for the strains of interest output file
    kmer_stroi.write("cluster\tstrain\tfeature_id\tcontig\tfeature_strand\tcontig_start\tcontig_end\tgene_start\tgene_end\tstrand\tk-mer\n")
    kmer_stroi.flush()

    return kmer_stroi


def create_hash_files(output, compress=False):
    if not compress:
        hash_pat = open(os.path.join(output, "hashes_to_patterns.tsv"), "w")
        kmer_hash = open(os.path.join(output, "kmers_to_hashes.tsv"), "w")
    else:
        hash_pat = gzip.open(os.path.join(output, "hashes_to_patterns.tsv.gz"),
                             "wt", compresslevel=9)
        kmer_hash = gzip.open(os.path.join(output, "kmers_to_hashes.tsv.gz"),
                              "wt", compresslevel=9)

    return hash_pat, kmer_hash


def create_faidx(file_name, uppercase=True, close=True):

    faidx = Fasta(file_name,
                  sequence_always_upper=uppercase,
                  rebuild=False)
    if close:
        faidx.close()
        return file_name
    else:
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


def iter_gene_clusters(panaroo, genome_data, up, down, down_start_codon, patfilt,
                       gene_list=None, raise_missing=False):
    # check if we have genome data for all
    # strains in the pangenome
    # if not give a warning
    all_strains = set(panaroo.columns)
    missing_strains = all_strains.difference(genome_data.keys())
    if len(missing_strains) > 0:
        logger.warning(f'There are {len(missing_strains)} strains present '
                       'in the pangenome table but not in the GFF directory')
        for strain in missing_strains:
            logger.debug(f'Missing GFF file: {strain}')
        if raise_missing:
            raise KeyError(f'Missing {len(missing_strains)} from the GFF directory')

    # go through each gene cluster
    all_ogs = panaroo.shape[0]
    for i, (idx, row) in enumerate(panaroo.iterrows()):
        if gene_list is not None and idx not in gene_list:
            logger.debug(f"Skipping {idx} ({i+1}/{all_ogs})")
            continue

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

        clusterpresab = np.zeros(len(strains), dtype = int)
        for strain in present:
            clusterpresab[sortstrain[strain]] = 1

        # print(clusterpresab)
        # cycle through all the strains that have the gene
        for strain, genes in row.dropna().items():

            strain = str(strain)

            if strain not in genome_data:
                continue

            gene_sequences[strain] = []
            # access the GFF/fasta data
            sequences_fname, features = genome_data[strain]
            sequences = create_faidx(sequences_fname, close=False)
            # be aware of paralogs
            for gene in genes.split(';'):
                # access a particular gene
                try:
                    feat = features[gene]
                except KeyError:
                    # e.g. refound genes are not in the GFF
                    logger.warning(f"Could not find gene {gene} from {idx} in {strain}")
                    if raise_missing:
                        raise KeyError(f"Could not find gene {gene} from {idx} in {strain}")
                    continue
                # double check if the chromosome is in the fasta file
                try:
                    sequences[feat.chromosome]
                except KeyError:
                    # e.g. refound genes are not in the GFF
                    logger.warning(f"Could not find chromosome {feat.chromosome} in {strain}")
                    if raise_missing:
                        raise KeyError(f"Could not find chromosome {feat.chromosome} in {strain}")
                    continue

                # corner case: upstream offset is over the contig's edge
                # adjust if so and save the true offset
                if feat.strand > 0 and feat.start-1-up < 0:
                    offset = feat.start - 1
                else:
                    offset = up
                # corner case: downstream offset is over the contig's edge
                # adjust if so and save the true offset
                if feat.strand < 0 and feat.start-1-down < 0:
                    offset_d = feat.start - 1
                else:
                    offset_d = down

                # access its gene sequence
                if not down_start_codon:
                    # down means relative to the gene's 3'
                    if feat.strand > 0:
                        seq = sequences[feat.chromosome][feat.start-1-offset:feat.end+offset_d]
                        seq_start = feat.start - offset
                        seq_end = feat.end + offset_d
                    else:
                        seq = -sequences[feat.chromosome][feat.start-1-offset_d:feat.end+offset]
                        seq_start = feat.start - offset_d
                        seq_end = feat.end + offset
                else:
                    # down means relative to the gene's 5'
                    if feat.strand > 0:
                        seq = sequences[feat.chromosome][feat.start-1-offset:feat.start+offset_d]
                        seq_start = feat.start - offset
                        seq_end = feat.start + offset_d
                    else:
                        seq = -sequences[feat.chromosome][feat.end-1-offset_d:feat.end+offset]
                        seq_start = feat.end - offset_d
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

            # close faidx handle
            sequences.close()

        # provide an empty list if the strain does not have the gene
        for strain in absent:
            gene_sequences[strain] = []

        yield gene_sequences, idx, clusterpresab


if __name__ == "__main__":
    pass
