# panfeed

panfeed makes use of panaroo's output files (gene_presence_absence.csv and gene_data.csv), which already cluster the genes of the input genomes into groups of close relatives.
panfeed then prepares the data as input for pyseer, while keeping the positional information for each k-mer that originates from a strain of interest (as specified by user), to facilitate the visualization and other downstream applications.

	python main.py  --stroi_in stroi.txt --stroi_out k-mers_for_strains_of_interest.csv --kmer_length 31 --presence_absence gene_presence_absence.csv --gene_data "" --start_inter 0 --end_inter 0 --canon True --specfilt True

# Documentation
To start with panfeed, simply copy the three scripts (ENTER SCRIPT NAMES) into the folder, containing the GFF files and the output files from panaroo. Then enter the command above (with your specific file names etc.).
The output comes in the form of 3 CSV files (+one error_log file that keeps track of refound genes), which contain information on the k-mers from strains of interest, hashed patterns to patterns and k-mers to hashed patterns. All files are tab delimited CSV files.

	python main.py  --stroi_in		#May be a simple textfile that contains the strain names for the strains of interest. K-mer positions for these strains are being recorded.
                                                                  
                	--stroi_out		#File name for the outputfile, which contains the logged k-mer positions. 
                
                	--kmer_length		#K-mer length that should be used by panfeed.
                
                	--presence_absence	#The presence_absence_gene_matrix as output by panaroo. Contains the cluster information for the different strains.
                                                                  
                	--gene_data		#This matrix contains additional cluster information, such as the refound genes and their sequence (but not their positional information).
                                                                  
                	--start_inter		#Includes the specified amount of nucleotides upstream of all genes (to include the promotor region for instance).
                
                	--end_inter		#Includes the specified amount of nucleotides downstream of all genes.
                
                	--canon			#Makes panfeed only consider canonical k-mers (lexicographically).
                
                	--specfilt		#Filters out k-mers with patterns that correspond to the original pattern of the entire cluster. If this pattern is significant for a certain trait, it is likely due to the presence or absence of the entire gene.

# Prerequisites:
`pyfaidx` (0.6.3.1)

`numpy` (1.20.3)

`pandas` (1.3.2)

# Citations
`pyseer`: Lees, John A., Galardini, M., et al. pyseer: a comprehensive tool for microbial pangenome-wide association studies. Bioinformatics 34:4310–4312 (2018). doi:10.1093/bioinformatics/bty539.

`panaroo`: Tonkin-Hill G, MacAlasdair N, Ruis C, Weimann A, Horesh G, Lees JA, Gladstone RA, Lo S, Beaudoin C, Floto RA, Frost SDW, Corander J, Bentley SD, Parkhill J. 2020. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol 21:180.
