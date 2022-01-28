# panfeed

panfeed makes use of panaroo's output file (gene_presence_absence.csv), which already clusters the genes of the input genomes into groups of close relatives.
panfeed then prepares the data as input for pyseer, while keeping the positional information for each k-mer that originates from a strain of interest (as specified by user), to facilitate the visualization and other downstream applications.

# Prerequisites:

`pyfaidx` (0.6.3.1)

`numpy` (1.20.3)

`pandas` (1.3.2)

# Getting started
Download the zipped folder and unpack it in the desired directory. You can then either install panfeed or use the panfeed-runner.py script to run it.

Installation:

	python .\setup.py install

Panfeed may then be called via the console, irrespective of the directory. e.g.:

	panfeed -g ./gffdirectory/ --targets ./target.txt -o ./panfeedoutput -p gene_presence_absence.csv


To use panfeed via the runner script, simply execute the script by typing:

	python ./panfeed-runner.py
	
(the path is of course dependent on your current position)



# Documentation
	-h
Displays the following commands and their description.

	-g --gff
Directory which contains all the GFF files. They must also contain the genome sequence and should be named the same as the panaroo .csv-header

	--targets
File indicating, for which samples the k-mer positions should be logged. a simple text file containing the sample names suffices.

	-o --output
Name of the output directory. If the directory already exists, an error will be raised.

	-k --kmer-length
Length of the k-mers. The default is 31 nt. 

	-p --presence-absence
Panaroo's gene_presence_absence.csv, which contains all the clustering information.

	--maf
Set the minor allele frequency threshold. K-mers that are represented in a percentage of genomes below this threshold are excluded. --maf must be below 0.5 (default = 0.1).

	--upstream
How many base pairs to include upstream of the actual gene sequence. (e.g. to include promoter regions)

	--downstream
How many base pairs to include downstream of the actual gene sequence. (e.g. to include promoter regions)

	--downstream-start-codon
If active, the --downstream parameter will be based on the start codon coordinates of the gene.

	--non_canonical
This option forces panfeed to also compute the non-canonical k-mers. By default, panfeed only considers the canonical (lexicographically smallest) k-mers.

	--no-filter
By default panfeed filters out all k-mers that have the same presence absence pattern as the gene cluster. This option may be activated to include these k-mers.

	--multiple-files
Activating this option will generate one set of output files for each gene cluster.

	-v
If this option is used, the verbosity of the displayed information is increased.

	--version
Displays the current version of panfeed.

The output of panfeed comes in the form of three TSV files that contain information on the position of logged k-mers (kmers.tsv), hashed patterns to patterns (hashes_to_patterns.tsv) and k-mers to hashed patterns (kmers_to_hashes.tsv). Additionally, panfeed will produce the .fasta and .fasta.fai (fasta index) files for the respective genomes/GFF files.

# Citations
`pyseer`: Lees, John A., Galardini, M., et al. pyseer: a comprehensive tool for microbial pangenome-wide association studies. Bioinformatics 34:4310–4312 (2018). doi:10.1093/bioinformatics/bty539.

`panaroo`: Tonkin-Hill G, MacAlasdair N, Ruis C, Weimann A, Horesh G, Lees JA, Gladstone RA, Lo S, Beaudoin C, Floto RA, Frost SDW, Corander J, Bentley SD, Parkhill J. 2020. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol 21:180.
