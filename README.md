# panfeed

[![Anaconda-Server Badge](https://anaconda.org/bioconda/panfeed/badges/version.svg)](https://anaconda.org/bioconda/panfeed)

`panfeed` is a k-mer streaming tool that works one gene cluster at a time.
Starting from a list of annotated genome assemblies in GFF3 format and a gene presence absence matrix (as produced by
[`roary`](https://github.com/sanger-pathogens/Roary), [`panaroo`](https://github.com/gtonkinhill/panaroo/)
and [`ggCaller`](https://github.com/samhorsfield96/ggCaller)),
`panfeed` generates a table with unique k-mer presence/absence patterns,
which can be used for genome-wide associations (GWAS) using tools
such as [`pyseer`](https://github.com/mgalardini/pyseer).
Mapping of associated patterns to gene clusters and base resolution mapping of k-mers can be then achieved with the other two outputs of `panfeed`.
Advantages of this approach over the generation of k-mers from a global de Bruijn
graph include a lower chance of encountering artifacts
due to repetitive regions and easier interpretations and visualization of results.

# Citation

Sommer, H., Djamalova, D., & Galardini, M. (2023). **Reduced ambiguity and improved interpretability of bacterial genome-wide associations using gene-cluster-centric k-mers.** Microbial Genomics, 9(11), 001129. [10.1099/mgen.0.001129](https://doi.org/10.1099/mgen.0.001129)

# Installation

`panfeed` can be installed using `pip`:

    python3 -m pip install panfeed

Or through `conda` (or `mamba` to speed things up):

    conda create -n panfeed -c bioconda panfeed

Alternatively, we provide a `conda` recipe to create an environment 
named `panfeed `. Download the
[environment file](https://raw.githubusercontent.com/microbial-pangenomes-lab/panfeed/main/environment.yml)
and then run:

    conda env create -f environment.yml
    conda activate panfeed

# Quick start guide

We reccommend a two-pass approach when using `panfeed`; the first pass generates the
presence/absence matrix for all k-mers across all gene clusters. After the association
analysis is completed, the `panfeed-get-clusters` command can be used to list the gene
clusters with k-mers passing the desired significance threshold, and a second pass of
`panfeed` can be run on those gene clusters alone to generate a base-level mapping of
all k-mers across samples for fine-mapping and visualization purposes.
The main advantage of the two-pass approach is a significant reduction in storage
requirements for the k-mer metadata file, and a slightly shorter computation time.

To run the first pass, prepare a folder with all GFF3 annotated assemblies files
(including the nucleotide sequences at the end of each file), with
file name in the format `SAMPLE.gff`. Each sample name should have a matching column
in the gene clusters presence/absence file, which must follow the same format as those
generated by `panaroo` (_i.e._ `gene_presence_absence.csv`).
Then run the following command, which will include 100 bases upstream and
downstream of each gene cluster:

    panfeed -g gffs -p gene_presence_absence.csv -o panfeed1 --upstream 100 --downstream 100 --compress --cores 4

This will create three files in the `panfeed1` directory:

* `kmers.tsv.gz`: k-mers metadata file (empty for this pass)
* `kmers_to_hashes.stv.gz`: file to match gene clusters, k-mer sequences and the hash for the respective presence/absence pattern
* `hashes_to_patterns.tsv.gz`: binary presence/absence matrix for all unique k-mer patterns (rows) across samples (columns)

The `hashes_to_patterns.tsv.gz` file can be used to run a GWAS analysis
with a tool such as `pyseer`, which will produce an output table
(_e.g._ `pyseer.tsv`) with association
statistics for each pattern passing the basic filtering thresholds. This file can
then be used to retrieve the gene clusters that encode k-mers passing the desired
significance threshold:

    panfeed-get-clusters -a pyseer.tsv -p panfeed1/kmers_to_hashes.stv.gz -t 1E-7 > gene_clusters.txt

The second pass of `panfeed` can be then run focusing on the "interesting"
gene clusters and generating k-mers positional information across all samples:

    ls gffs/ | sed 's/.gff//g' > samples.txt
    panfeed -g gffs -p gene_presence_absence.csv -o panfeed2 --targets samples.txt --genes gene_clusters.txt --upstream 100 --downstream 100 --compress --cores 4

This time the `kmers.tsv.gz` file will contain absolute and relative (to start codon) positional information
for each k-mer in all samples. This file can be then merged with the association results so that the
association statistics are paired with each k-mer and their position across samples:

    panfeed-get-kmers -a pyseer.tsv -p kmers_to_hashes.tsv.gz -k kmers.tsv.gz | gzip > annotated_kmers.tsv.gz

Association results can be then visualized for each gene cluster with the following command,
which requires the phenotype file used for the association with `pyseer` (`data.tsv`
in the example command below):

    panfeed-plot -k annotated_kmers.tsv.gz -p data.tsv

This command will generate three figures for each gene cluster:

* `significance_CLUSTER.png`: k-mers colors are proportional to their significance level, think one Manhatten plot for each sample
* `sequence_CLUSTER.png`: k-mers are colored based on their nucleotide sequence, effectively generating a pseudo-alignment
* `hybrid_CLUSTER.png`: a combination of the two previous figures; k-mers color is based on their sequence, opacity is proportional to their significance level

Additionally, a file called `sequence_legend.png` is created to indicate which color is associated to which nucleotide.

# Additional information

If your GFF files do not contain the nucleotide sequences, you can provide them to `panfeed`
as a separate argument, using `-f fastas`. The `fastas` folder should contain one file per sample
with the name format `SAMPLE.fasta` or `SAMPLE.fna`.

If you want the k-mers presence/absence patterns to encode differently the information on
whether a gene cluster is missing from a sample, use the `--consider-missing` argument.
By default a missing gene cluster is encoded as `0`, same as a missing k-mer.

The visualization command has many arguments to fully customize the resulting plots;
among them:

* `--phenotype-column PHENOTYPE`, will sort the plots by the provided phenotype value (descending order)
* `--start -50 --stop 100 --sample 0.1`, will restrict the plot to 10% of samples and to the -50 to +100 region relative to the start codon
* adding `--nucleotides` to the above command will add the nucleotide letters to each plot

# Prerequisites:

The following packages and version have been used to develop and test `panfeed`

* `pyfaidx` (0.6.3.1)
* `numpy` (1.20.3)
* `pandas` (1.3.2)
* `matplotlib` (3.5.2)
* `seaborn` (0.11.2)

