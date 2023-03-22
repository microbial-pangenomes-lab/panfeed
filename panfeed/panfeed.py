#!/usr/bin/env python

import os
import logging
import numpy as np
import hashlib
import binascii
from copy import deepcopy
from io import StringIO

from .input import create_hash_files, create_kmer_stroi

logger = logging.getLogger('panfeed.panfeed')


def init_presabs_vector(n_strains, clusterpresab, missing_nan=False):
    v = np.zeros(n_strains, dtype=np.float64)
    if missing_nan:
        v[clusterpresab == 0] = np.nan
    return v


def cluster_cutter(cluster_gen, klength, stroi,
                   multiple_files, canon, consider_missing_cluster,
                   output, compress=False):
#iterates through genes present in the cluster
#shreds gene sequence into k-mers and adds them to the cluster dictionary
#if a gene belongs to a strain of interest, additional info on the k-mer is saved
    kmer_stroi = None

    cluster, idx, clusterpresab = cluster_gen
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
        kmer_stroi = create_kmer_stroi(path, compress)

    cluster_dict = {}

    sortstrain = {x: i
                    for i, x in enumerate(sorted(cluster.keys()))}
    n_strains = len(sortstrain)

    base_vector = init_presabs_vector(n_strains, clusterpresab,
                                      consider_missing_cluster)

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
                        used_strand = 1
                    else:
                        canonseq = revspecseq
                        used_strand = -1 

                    cluster_dict[canonseq] = cluster_dict.get(canonseq,
                                                              deepcopy(base_vector))
                    cluster_dict[canonseq][sortstrain[strain]] = 1
                else:
                    used_strand = strand
                    cluster_dict[specseq] = cluster_dict.get(specseq,
                                                             deepcopy(base_vector))
                    cluster_dict[specseq][sortstrain[strain]] = 1

                    cluster_dict[revspecseq] = cluster_dict.get(revspecseq,
                                                                deepcopy(base_vector))
                    cluster_dict[revspecseq][sortstrain[strain]] = 1

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
                        memchunk.write(f"{idx}\t{strain}\t{gene_id}\t{contig}\t{strand}\t{truestart}\t{trueend}\t{genestart}\t{geneend}\t{used_strand}\t{canonseq}\n")
                    else:
                        memchunk.write(f"{idx}\t{strain}\t{gene_id}\t{contig}\t{strand}\t{truestart}\t{trueend}\t{genestart}\t{geneend}\t{used_strand}\t{specseq}\n")
                        memchunk.write(f"{idx}\t{strain}\t{gene_id}\t{contig}\t{strand}\t{truestart}\t{trueend}\t{genestart}\t{geneend}\t{-used_strand}\t{revspecseq}\n")

    if multiple_files:
        kmer_stroi.write(memchunk.getvalue())
        return idx, cluster_dict, clusterpresab, None
    else:
        return idx, cluster_dict, clusterpresab, str(memchunk.getvalue())
 

def write_headers(hash_pat, kmer_hash, genepres):
    memchonkheader1 = StringIO()
    memchonkheader2 = StringIO()
    
    memchonkheader1.write("hashed_pattern")
    for strain in sorted(genepres.columns):
        memchonkheader1.write(f"\t{strain}")
    memchonkheader1.write("\n")
    hash_pat.write(memchonkheader1.getvalue())
    hash_pat.flush()

    memchonkheader2.write("cluster\tk-mer\thashed_pattern\n")
    kmer_hash.write(memchonkheader2.getvalue())
    kmer_hash.flush()


def pattern_hasher(cluster_dict_iter, kmer_stroi, hash_pat, kmer_hash,
                   genepres, patfilt, maf, output, patterns=None,
                   consider_missing_cluster=False,
                   compress=False):
    #iterates through the cluster dictionary output by cluster_cutter()
    #outputs the hashed patterns, patterns and k-mers to the output files
    #two files are created: hashed k-mer patterns to presence/absence patterns
    #                       k-mer sequences to hashed k-mer patterns
    
    # set flag for multiple files
    multiple_files = False
    if hash_pat is None or kmer_hash is None:
        multiple_files = True
    
    # keep track of already observed patterns
    # might have a big memory footprint
    # TODO: check for memory footprint
    if patterns is None:
        patterns = set()

    for idx, cluster_dict, clusterpresab, memchunk in cluster_dict_iter:
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
            hash_pat, kmer_hash = create_hash_files(path, compress)
            # must "flush" the already observed patterns file
            patterns = set()

            write_headers(hash_pat, kmer_hash, genepres)
        else:
            if memchunk is not None:
                kmer_stroi.write(memchunk)

        memchunkhash_pat = StringIO()
        memchunkkmer_hash = StringIO()
       
        hashed = hashlib.md5(clusterpresab.view(np.uint8))
        khash = binascii.b2a_base64(hashed.digest()).decode()[:24]
        memchunkkmer_hash.write(f"{idx}\t\t{khash}\n")
        
        if khash not in patterns:
            patterns.add(khash)
            if not consider_missing_cluster:
                patterntup = "\t".join(map(str, clusterpresab.astype(np.uint8)))
            else:
                patterntup = "\t".join(['' if np.isnan(x)
                                        else str(int(x))
                                        for x in clusterpresab])
            memchunkhash_pat.write(f"{khash}\t{patterntup}\n")
        
        for kmer in cluster_dict:
            if not consider_missing_cluster:
                af = cluster_dict[kmer].sum() / cluster_dict[kmer].shape[0]
            else:
                # be aware of missing values
                k_presab = cluster_dict[kmer]
                k_presab = k_presab[~np.isnan(k_presab)]
                af = k_presab.sum() / k_presab.shape[0]
            if af >= 0.5:
                af = 1-af
            if af < maf:
                continue
            
            if patfilt == True:
                if tuple(cluster_dict[kmer]) == tuple(clusterpresab):
                    continue

            hashed = hashlib.md5(cluster_dict[kmer].view(np.uint8))
            khash = binascii.b2a_base64(hashed.digest()).decode()[:24]
            memchunkkmer_hash.write(f"{idx}\t{kmer}\t{khash}\n")
            
            if khash in patterns:
                continue
            patterns.add(khash)

            if not multiple_files and not len(patterns) % 1000:
                logger.debug(f"Observed patterns so far: {len(patterns)}")

            if not consider_missing_cluster:
                patterntup = "\t".join(map(str, cluster_dict[kmer].astype(np.uint8)))
            else: 
                patterntup = "\t".join(['' if np.isnan(x)
                                        else str(int(x))
                                        for x in cluster_dict[kmer]])
            memchunkhash_pat.write(f"{khash}\t{patterntup}\n")

        hash_pat.write(memchunkhash_pat.getvalue())
        kmer_hash.write(memchunkkmer_hash.getvalue())

        if multiple_files:
            memchunkhash_pat = None
            memchunkkmer_hash = None

    hash_pat.flush()
    kmer_hash.flush()

    return patterns

if __name__ == "__main__":
    pass
