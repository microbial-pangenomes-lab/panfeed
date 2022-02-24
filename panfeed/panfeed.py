#!/usr/bin/env python

import os
import logging
import numpy as np
import hashlib
import binascii
from io import StringIO

from .input import create_hash_files, create_kmer_stroi

logger = logging.getLogger('panfeed.panfeed')

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
                            cluster_dict[f"{canonseq}"] = np.zeros(n_strains, dtype = int)
                        cluster_dict[f"{canonseq}"][sortstrain[strain]] = 1
                    else:
                        used_strand = strand
                        if specseq not in cluster_dict:
                            cluster_dict[f"{specseq}"] = np.zeros(n_strains, dtype = int)
                        cluster_dict[f"{specseq}"][sortstrain[strain]] = 1
                        
                        if revspecseq not in cluster_dict:
                            cluster_dict[f"{revspecseq}"] = np.zeros(n_strains, dtype = int)
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

        if multiple_files:
            kmer_stroi.write(memchunk.getvalue())
            memchunk = None
         
        yield idx, cluster_dict, clusterpresab, memchunk
 

def write_headers(hash_pat, kmer_hash, genepres):
    memchonkheader1 = StringIO()
    memchonkheader2 = StringIO()
    
    memchonkheader1.write("hashed_pattern")
    for strain in sorted(genepres.columns):
        memchonkheader1.write(f"\t{strain}")
    memchonkheader1.write("\n")
    hash_pat.write(memchonkheader1.getvalue())

    memchonkheader2.write("k-mer\thashed_pattern\n")
    kmer_hash.write(memchonkheader2.getvalue())


def pattern_hasher(cluster_dict_iter, kmer_stroi, hash_pat, kmer_hash, genepres, patfilt, maf, output):
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
            hash_pat, kmer_hash = create_hash_files(path)
            # must "flush" the already observed patterns file
            patterns = set()

            write_headers(hash_pat, kmer_hash, genepres)
        else:
            if memchunk is not None:
                kmer_stroi.write(memchunk.getvalue())

        memchunkhash_pat = StringIO()
        memchunkkmer_hash = StringIO()
        
        pattern = clusterpresab.view(np.uint8)
        hashed = hashlib.md5(pattern)
        khash = binascii.b2a_base64(hashed.digest()).decode()[:24]
        memchunkkmer_hash.write(f"{idx}\t{khash}\n")
        
        if khash not in patterns:
            patterns.add(khash)
            patterntup = "\t".join(map(str, clusterpresab))
            memchunkhash_pat.write(f"{khash}\t{patterntup}\n")
        
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

        if multiple_files:
            hash_pat.write(memchunkhash_pat.getvalue())
            kmer_hash.write(memchunkkmer_hash.getvalue())
            memchunkhash_pat = None
            memchunkkmer_hash = None


if __name__ == "__main__":
    pass
