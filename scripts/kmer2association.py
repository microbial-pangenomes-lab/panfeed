from collections import defaultdict, namedtuple
import shutil
import os
from io import StringIO
import argparse

def get_options():
    description = """script that allows the correlation of GWAS associations 
                    and meta data provided by panfeed"""
    
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument("-p", "--presence_absence",
                        help="CSV file containing the presence/absence matrix as output by panaroo")

    parser.add_argument("-g", "--gffs",
                        help="folder containing the GFF files")

    parser.add_argument("-asso", "--associations_file",
                        help="TSV file containing hashes and their significance (pyseer output)")
    
    parser.add_argument("-kmers", "--kmer_file",
                        help="TSV file containing kmers and related info")
    
    parser.add_argument("-kth", "--kmers_to_hashes",
                        help="TSV file relating kmers and their hashes")
    
    parser.add_argument("-htp", "--hashes_to_patterns",
                        help="TSV file relating hashes and patterns")
    
    parser.add_argument("-out", "--out_path", 
                        help="path where output is thrown")
    
    parser.add_argument("-alpha", "--alpha", type = float,
                        default=0.05,
                        help="family  wise error-rate")
    
    return parser.parse_args()
   

Feature = namedtuple('Feature', ['id', 
                                 'chromosome',
                                 'start', 
                                 'end', 
                                 'strand'])

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
                continue
    
    return features

def det_thresh(htpfile, alpha):
    #count unique patterns to calculate association significance threshold 
    with open(htpfile, "r") as htp:
        linecount = 0
        next(htp)
        for line in htp:
            linecount += 1
        
        #remove count for trailing "\n"
        linecount -= 1
    
    threshold = alpha / float(linecount)
    
    return threshold

def sig_sort(assofile, threshold, outp):
    #filter out insignificant associations according to previously determined threshold
    sigdict = {}
    buffr = StringIO()
    
    with open(assofile, "r") as asfi:
        next(asfi)

        for line in asfi:
            splits = list(line.rstrip("\n").split("\t"))

            #apply threshold to lrt-pvalue of associations file
            if float(splits[3]) <= threshold:
                sigdict[splits[0]] = splits[1:9]
    
    with open(assofile, "r") as asfi:
        header = str(asfi.readline())

    buffr.write(header)
    filtfile = open(f"{outp}filtered.tsv", "w")

    #write all the filtered associations to filtered.tsv
    for key in sigdict:
        pyseerdata = '\t'.join(sigdict[key])
        buffr.write(f"{key}\t{pyseerdata}\n")

    filtfile.write(buffr.getvalue())
    filtfile.flush()
    filtfile.close()

    return sigdict
    
def kmer_corr(sigdict, hash_kmer, outp):
    #link k-mers to their respective hashes (/patterns)
    kmerdict = {}
    clusterset = set()

    with open(hash_kmer, "r") as hakm:
        next(hakm)
        
        for line in hakm:
            splits = line.rstrip("\n").split("\t")
            
            if splits[2] in sigdict:
                kmerdict[f"{splits[0]}:{splits[1]}"] = splits[2] #explanation: splits[0] = cluster, splits[1] = k-mer, splits[2] = hashed pattern
                clusterset.add(splits[0])

    with open(f"{outp}sig_clusters.txt", "w") as sicl:
        for clust in clusterset:
            sicl.write(f"{clust}\n")

    #kmerdict contains {"cluster:kmer": hashedpattern, "cluster:kmer": hashedpattern}
    return kmerdict

def kmer_data_picker(sigdict, kmerdict, kmerfile, outp):
    #write all significant (and logged) k-mers to file (with all linked information)
    bigbuffer = StringIO()
    duplset = set()

    #header for sig_kmer_data.tsv
    bigbuffer.write("k-mer\thashed\tcluster\tstrain\tfeature_id\tcontig\tfeature_strand\tcontig_start\tcontig_end\tgene_start\tgene_end\tstrand\taf\tfilter-pvalue\tlrt-pvalue\tbeta\tbeta-std-err\tvariant_h2\tnotes\n")
    with open(kmerfile, "r") as kfile:

        for line in kfile: #iterate through kmers.tsv
            splits = line.rstrip("\n").split("\t")
            kmer = splits[10]
            cluster = splits[0]

            if f"{cluster}:{kmer}" in kmerdict: #look for significant cluster:k-mer combinations in kmers.tsv
                kmerdictkey = f"{cluster}:{kmer}"
                panfeeddata = '\t'.join(splits[0:10])
                beddata = (splits[3], splits[5], splits[6], splits[9]) #information that is included in bed files (contig, start_of_sequence, end_of_sequence, strand)

                if beddata not in duplset: #write to buffer, only if the unique sequence position has not yet been found in another cluster
                    duplset.add(beddata)
                    pyseerdata = '\t'.join(sigdict[kmerdict[kmerdictkey]])
                    bigbuffer.write(f"{kmer}\t{kmerdict[kmerdictkey]}\t{panfeeddata}\t{pyseerdata}\n")

                else:
                    continue

    return bigbuffer

def cluster_data_fetcher(kmerdict, sigdict, gffp, presab, bigbuffer, outp):
    #finds clusters in which the presence/absence of the entire gene is significant for the phenotype
    clusterlist = []

    for clustmer in kmerdict:
        if clustmer.split(":")[1] == "":
            clusterlist.append(clustmer.split(":")[0])

    clusterdict = defaultdict(list)

    with open(f"{presab}", "r") as paf:
        header = paf.readline().split(",")[3:]

        for row in paf:
            splitrow = row.rstrip().split(",")

            for cluster in clusterlist: #iterates through the list of significant clusters
                if f"{splitrow[0]}" == cluster:
                    genenames = splitrow[3:] #all gene names within a cluster
            
                    for genename in genenames:
                        if genename != "nan" and genename != "" and "refound" not in genename:
                            if ";" in genename: #looks for strains with multiple gene paralogs in one cluster
                                paralogs = genename.split(";")

                                for paralog in paralogs:
                                    clusterdict[splitrow[0]].append([header[genenames.index(genename)], paralog])

                            else:    
                                clusterdict[splitrow[0]].append([header[genenames.index(genename)], genename])
                        else:
                            pass

    #clusterdict contains {clustername: [[strainname1, genename], [strainname2, genename]], ...}
    for cluster in clusterdict:

        for item in clusterdict[cluster]: #iterates through genes in the significant clusters and fetches the positional information from the GFF files
            cleanstrainname = item[0].rstrip("\n")
            dictgff = parse_gff(f"{gffp}{cleanstrainname}.gff")
            hashed = str(kmerdict[cluster+":"])
            pyseerdata = '\t'.join(sigdict[hashed])
            bigbuffer.write(f"\t{hashed}\t{cluster}\t{item[0].rstrip()}\t{item[1]}\t{dictgff[item[1]].chromosome}\t{dictgff[item[1]].strand}\t{dictgff[item[1]].start}\t{dictgff[item[1]].end}\t\t\t{dictgff[item[1]].strand}\t{pyseerdata}\n")

    with open(f"{outp}sig_kmer_data.tsv", "w") as sigkmer:
        sigkmer.write(bigbuffer.getvalue())
        sigkmer.flush()

if __name__ == "__main__":
    #get arguments
    args = get_options()
    
    outp = args.out_path
    kth = args.kmers_to_hashes
    kmerfile = args.kmer_file
    assofile = args.associations_file
    htp = args.hashes_to_patterns
    alpha = args.alpha
    gffs = args.gffs
    presab = args.presence_absence
    
    #check for output directory
    if not os.path.isdir(outp):
    	os.mkdir(out)

    #determine threshold for association significance (alpha / number of unique patterns)
    thresh = det_thresh(htp, alpha)
    
    #find all significant entries in association.tsv based on previous threshold and save significant patterns as filtered.tsv
    sigdict = sig_sort(assofile, thresh, outp)
    
    #link k-mers to their corresponting patterns/hashes
    kmerdict = kmer_corr(sigdict, kth, outp)
    
    #links logged k-mer data to the relevant k-mers and outputs sig_kmer_data.tsv 
    bigbuffer = kmer_data_picker(sigdict, kmerdict, kmerfile, outp)

    #fetches the sequence information for the significant clusters 
    #(where the absence/presence of the entire gene is significant) 
    #and writes them to file
    cluster_data_fetcher(kmerdict, sigdict, gffs, presab, bigbuffer, outp)


    
    
    
    
