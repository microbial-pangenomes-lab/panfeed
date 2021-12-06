import os, sys
import pandas as pd
from pyfaidx import Fasta, Faidx
from marcos_modules import parse_gff, create_faidx
import numpy as np
import hashlib
import binascii
import shutil
from io import StringIO
from re import sub

def what_are_my_inputfiles():
    filelist = []
    
    for file in os.listdir():
        
        if file.endswith(".gff" or "gff"):
            
            filelist.append(file.split(".")[0])
            
    if filelist == []:
        
        sys.stderr.write("no GFF files found in working directory :(")
        
        sys.exit(1)
        
    return filelist


def prep_data_n_fasta(filelist):#TODO edit to get input from user?
    
    data = {}
    
    for genome in filelist:
    
        txtfile = open(f"{genome}.fasta", "w")
    
        txtfile.write(open(f"{genome}.gff", "r").read().split("##FASTA")[1])
    
        txtfile.close()
    
        sequences = create_faidx(f'{genome}.fasta')
    
        features = parse_gff(f'{genome}.gff')
    
        data[genome] = (sequences, features)
    
    return data

def set_input_output(stroi_in, stroi_out, presence_absence, gene_data):
    #determines the names of the input, output files 
    
    if gene_data == "":
        
        pass
    
    else:
        
        gene_data = pd.read_csv(f"{gene_data}", sep=",")# load the gene_data_table that is needed to locate the sequence of refound genes
 
    genepres = pd.read_csv(f"{presence_absence}", sep=",", index_col=0).drop(columns=['Non-unique Gene name', 'Annotation'])# load the gene_presence_absence matrix which will be used to ascribe strain/genes to clusters
    if stroi_in != "":
        
        stroi = open(f"{stroi_in}", "r").read()
        
    else:
        
        stroi = ""
    
    if os.path.exists("./output"):
        
        shutil.rmtree("./output")
        
        os.mkdir("./output")
        
    else:
        
        os.mkdir("./output")
    
    error_log = open("./output/unfound_genes.out", "a")

    kmer_stroi = open(f"./output/{stroi_out}", "a")
    #creates the header for the strains of interest output file
    kmer_stroi.write("cluster\tstrain\tcontig\tcontig_start\tcontig_end\tgene_start\tgene_end\tstrand\tk-mer\n")
    
    hash_pat = open("./output/hashes_to_patterns.out", "a")

    kmer_hash = open("./output/kmers_to_hashes.out", "a")

    return stroi, error_log, kmer_stroi, hash_pat, kmer_hash, gene_data, genepres




def cluster_cutter(cluster_gen, klength, stroi, kmer_stroi):
#iterates through genes present in the cluster
#shreds gene sequence into k-mers and adds them to the cluster dictionary
#if a gene belongs to a strain of interest, additional info on the k-mer is saved

    for cluster, idx in cluster_gen:#only needed if we do not multithread
        
        memchunk = StringIO()
        
        cluster_dict = {}
    
        sortstrain = {x: i
                      for i, x in enumerate(sorted(cluster.keys()))}
        n_strains = len(sortstrain)

        for strain in cluster.items():
        
            for seq in strain[1]:
        
                seqrevcomp = seq.reverse.complement
                # make the sequence a string to hopefully save some time
                seq1 = str(seq)
                seqrevcomp1 = str(seqrevcomp)

                #determines the number of k-mers in the sequence
                num_kmer = len(seq) - klength + 1
                
                seqlen = len(seq1)
                
                orient = strain[1][0].orientation
                
                contig = strain[1][0].name
                
                contigstart = int(seq.start)
                
                if strain[0] in stroi:
                    #if a gene belongs to a strain of interest the strain,contig,start,stop,strand,sequence
                    #of its k-mers are being written to kmer_stroi
                    for pos in range(num_kmer):
                    
                        specseq = str(seq1[pos:pos + klength])
                
                        specseqrev = str(seqrevcomp1[pos:pos + klength])
                
                        if specseq not in cluster_dict:
                            # init with an "empty" array
                            cluster_dict[f"{specseq}"] = np.zeros(n_strains).astype(int)
                    
                        # flip the "0" to "1"
                        cluster_dict[f"{specseq}"][sortstrain[strain[0]]] = 1
                
                        truestart = contigstart + pos
                                
                        trueend = truestart + klength
                
                        memchunk.write(f"{idx}\t{strain[0]}\t{contig}\t{truestart}\t{trueend}\t{pos}\t{pos+klength}\t{orient}\t{specseq}\n")
                
                        if specseqrev not in cluster_dict:
                            # init with an "empty" array
                            cluster_dict[f"{specseqrev}"] = np.zeros(n_strains).astype(int)
                    
                        # flip the "0" to "1"
                        cluster_dict[f"{specseqrev}"][sortstrain[strain[0]]] = 1
                        
                        #the start of the reverse complement is being saved 
                        #as the start of the reverse complement + the position of the iterator
                        #the end is the startposition + k-mer length
                        revstart = int(seqrevcomp.start) - pos
                                
                        revend = revstart - klength
                
                        #the strand of the reverse is always being saved as the opposit of the primary strand
                        orientrev = 1
                
                        if orient > 0:
                    
                            orientrev = -1
                    
                        memchunk.write(f"{idx}\t{strain[0]}\t{contig}\t{revstart}\t{revend}\t{seqlen-pos}\t{seqlen-pos-klength}\t{orientrev}\t{specseqrev}\n")
                    
                else:
                #does the same as the if-statement block on the same indent,
                #without logging the k-mer infos, as the strain is not a strain of interest
                    for pos in range(num_kmer):
                    
                        specseq = str(seq1[pos:pos + klength])
                
                        specseqrev = str(seqrevcomp1[pos:pos + klength])
        
                        if specseq not in cluster_dict:
                            # init with an "empty" array
                            cluster_dict[f"{specseq}"] = np.zeros(n_strains).astype(int)
                    
                        # flip the "0" to "1"
                        cluster_dict[f"{specseq}"][sortstrain[strain[0]]] = 1
                    
                        if specseqrev not in cluster_dict:
                            # init with an "empty" array
                            cluster_dict[f"{specseqrev}"] = np.zeros(n_strains).astype(int)
                    
                        # flip the "0" to "1"
                        cluster_dict[f"{specseqrev}"][sortstrain[strain[0]]] = 1
                        
        kmer_stroi.write(memchunk.getvalue())
        
        yield cluster_dict
        

def pattern_hasher(cluster_dict_iter, hash_pat, kmer_hash, genepres):
    #iterates through the cluster dictionary output by cluster_cutter()
    #outputs the hashed patterns, patterns and k-mers to the output files
    #two files are created: hashed k-mer patterns to presence/absence patterns
    #                       k-mer sequences to hashed k-mer patterns
    memchonkheader1 = StringIO()
    
    memchonkheader2 = StringIO()
    
    memchonkheader1.write("hashed_pattern")
    
    for strain in sorted(genepres.columns):
        
        memchonkheader1.write(f"\t{strain}")
    
    memchonkheader1.write("\n")
    
    hash_pat.write(memchonkheader1.getvalue())
    
    memchonkheader2.write("k-mer\thashed_pattern\n")
    
    kmer_hash.write(memchonkheader2.getvalue())
    
    for cluster_dict in cluster_dict_iter:
        
        memchunkhash_pat = StringIO()
        
        memchunkkmer_hash = StringIO()
        
        for kmer in cluster_dict:
        
            pattern = cluster_dict[kmer].view(np.uint8)
        
            hashed = hashlib.md5(pattern)
        
            khash = binascii.b2a_base64(hashed.digest()).decode()[:24]
        
            patterntup = "\t".join(map(str, cluster_dict[kmer]))
            
            memchunkhash_pat.write(f"{khash}\t{patterntup}\n")
            
            memchunkkmer_hash.write(f"{kmer}\t{khash}\n")
            
        hash_pat.write(memchunkhash_pat.getvalue())
        
        kmer_hash.write(memchunkkmer_hash.getvalue())






