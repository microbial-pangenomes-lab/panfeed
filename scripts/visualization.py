import os
from io import StringIO as sio
from statistics import mean
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import seaborn as sns
import copy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
from pyfaidx import Fasta
from math import isnan
from functools import partial
from matplotlib import colors
from matplotlib import rcParams
import matplotlib.patches as mpatches


############TO DO LIST:
		####replace unnecessary for loops with comprehensions (exchanging letters with numbers for seqdict)
		####look at end of script and remove arguments if not needed
		####


def get_options():
	description = """convenience script for creating plots 
	displaying the distribution of significant 
	k-mers in specific gene clusters"""
	
	parser = argparse.ArgumentParser(description=description)
	
	parser.add_argument("-phe", "--phenotypes",
						help="TSV file containing strains and their respective phenotype")

	parser.add_argument("--binary_phenotype", 
						action="store_true",
						default=False,
						help="""sets binary phenotype to True. 
								continuous phenotype is assumed by default.""")
	
	parser.add_argument("-sc", "--strain_column", type = int,
						help="column in phenotype.tsv containing the strain_ID")

	parser.add_argument("-pc", "--pheno_column", type = int,
						help="column in phenotype.tsv containing the phenotype")

	parser.add_argument("-pan", "--panaroo_matrix",
						help="panaroo output.csv (usually gene_presence_absence.csv)")

	parser.add_argument("-sikm", "--significant_kmers",
						help="CSV containing significant k-mers and associated data")

	parser.add_argument("-gff", "--gff_path",
						help="path to gff files")

	parser.add_argument("-htp", "--hashes2patterns",
						help="""CSV file containing patterns and
								their hashes, as output by panfeed""")

	parser.add_argument("-cid", "--cluster_ID",
						help="""Name of cluster for which to create the plot""")

	parser.add_argument("-alpha", "--family_wise_error_rate", type = float,
						default = 0.05,
						help="Alpha value that was used as pyseer input")

	parser.add_argument("-ud", "--upstream_downstream", nargs=2, type = int,
						default=[0,0],
						help="upstream and downstream region to include")

	parser.add_argument("-zoom", "--close_up_region", nargs=2, type = int,
						default=None,
						help="""default = None
								region to zoom into.
								must be None or two integers.
								positions are relative to the first position in the normal plot.""")

	parser.add_argument("--strain_marker", 
						action="store_true",
						default=False,
						help="""when active, indicates separation
								between positive and negative phenotypes.
								only use with binary phenotype!""")
	
	parser.add_argument("-outpath", "--output_path",
						default="./",
						help="path for output")
						
	parser.add_argument("-shift", "--shift_pattern", nargs=2, type = str,
						default=None,
						help="upstream and downstream region to include")

	parser.add_argument("-max", "--max_scalar", type = float,
						default=10.0,
						help="maximum scalar -log10 p-value to be displayed in the plot")

	parser.add_argument("--hybrid_heatmap", 
						action="store_true",
						default=False,
						help="""when active, nucleotides will be displayed as colors 
								and -log10 p-values as opacity.""")

	parser.add_argument("--nucleotides_only", 
						action="store_true",
						default=False,
						help="""when active, only the nucleotide sequence will be displayed
								in the resulting heatmap.""")

	parser.add_argument("--display_nucleotides", 
						action="store_true",
						default=False,
						help="""when active, nucleotides will be displayed in the plot.
								only recommended for small subsets of strains!""")

	parser.add_argument("-stroi", "--strains_of_interest",
						default=None,
						help="""file containing strains of interest. 
								must contain one strain ID per line.""")

	return parser.parse_args()

def calc_scalar(x, maxval, thresh):

	m = (1-0.2)/(maxval-thresh)

	if isnan(x) or type(x) is not np.float64 or x == 0. or x < threshold or x is []:
		y = 0.1

	elif x > maxval:
		y = 0.2 + m *(maxval - threshold)

	else:
		y = 0.2 + m *(x - threshold)

	return y

if __name__ == "__main__":

	args = get_options()

	phenofile = args.phenotypes
	phenotype = args.binary_phenotype
	straincol = args.strain_column
	phenocol = args.pheno_column
	panaroomat = args.panaroo_matrix
	sigkmers = args.significant_kmers
	gffpath = args.gff_path
	htpfile = args.hashes2patterns
	outpath = args.output_path
	clusterid = args.cluster_ID
	alpha = args.family_wise_error_rate
	upstream, downstream = args.upstream_downstream
	boundaries = args.close_up_region
	strainmarker = args.strain_marker
	shift = args.shift_pattern
	maxscalar = args.max_scalar
	hybrid = args.hybrid_heatmap
	nucleotides = args.nucleotides_only
	display_nucs = args.display_nucleotides
	stroi = args.strains_of_interest

	###check if upstream was assigned properly
	if upstream > 0:
		print("ERROR: upstream must be negative integer or 0")
		sys.exit()

	###check if hybrid_heatmap and nucleotides only are active at the same time
	if hybrid and nucleotides:
		print("ERROR: only one of --hybrid_heatmap and --nucleotides_only can be active at a time")
		sys.exit()

	###check if strainmarker is only active with binary_phenotype
	if not phenotype and strainmarker:
		print("ERROR: strainmarker can only be used with --binary_phenotype!")
		sys.exit()

	base2int = {'A': 0,
				'G': 1,
				'T': 2,
				'C': 3}

	float2legend = {"0.0": "A",
				"1.0": 'G',
				"2.0": 'T',
				"3.0": 'C',
				"nan": "end of contig"}

	#creating smaller versions of the panaroo presence/absence matrix and the phenotype file,
	#containing only the user specified strains.
	if stroi is not None:
		stroiset = set()
		with open(stroi, "r") as stroi:
			for line in stroi:
				stroiset.add(line.rstrip())

		with open(panaroomat, "r") as pama:
			delset = set(pama.readline().rstrip("\n").split(",")[3:])-stroiset

		newpresab = pd.read_csv(panaroomat, header=0, index_col=0, dtype = str)
		newpresab = newpresab.drop(columns=delset)
		newpresab.to_csv("temp_presab.csv")

		newphe = sio()
		with open(phenofile, "r") as phfi:
			newphe.write(phfi.readline())
			for line in phfi:
				if line.split("\t")[straincol] in stroiset:
					newphe.write(line)

		with open("temp_pheno.tsv", "w") as nfile:
			nfile.write(newphe.getvalue())

		panaroomat = "temp_presab.csv"
		phenofile = "temp_pheno.tsv"

	# iterating through presence/absence matrix until cluster of interest is found
	with open(panaroomat, "r") as pama: 
		strains = pama.readline().rstrip("\n").split(",")[3:]

		for line in pama:
			splitline = line.rstrip("\n").split(",")
			ID = splitline[0]

			if ID == clusterid:
				genenames = splitline[3:]
				break
	
	if phenotype:
		positivelist = []
		negativelist = []
		with open(phenofile, "r") as phefi:
			phefi.readline()
			for line in phefi:
				splili = line.rstrip().split("\t")
				strain = splili[straincol]
				phenot = int(splili[phenocol])

				if phenot == 1:
					positivelist.append(strain)
				elif phenot == 0:
					negativelist.append(strain)
				else:
					print("only binary phentypes allowed!")
					sys.exit()

	else:
		orderdict = {}
		with open(phenofile, "r") as phefi:
			phefi.readline()

			for line in phefi:
				splili = line.rstrip().split("\t")
				orderdict[splili[straincol]] = float(splili[phenocol])

	with open(sigkmers, "r") as sikm:
		sikm.readline()
		bestP = -np.log10(float(sikm.readline().rstrip("\n").split("\t")[14]))
		for line in sikm: # iterating through the significant k-mers file, trying to find the lowest lrt-p-value
			splitline = line.rstrip("\n").split("\t")
			logval = -np.log10(float(splitline[14]))
			if bestP < logval:
				bestP = logval
	
	seqdict = {}
	absdict = {}
	shiftdict = {}
	longestgene = 0
	paralogs = False
	ups = None
	downs = None
	
	for genename in range(len(genenames)): # iterating through the genenames in the cluster of interest
		strain = strains[genename]

		if ";" in genenames[genename]: # checking if paralogs are present or if a strain has no gene within the cluster
			names = genenames[genename].split(";")
			paralogs = True

		elif genenames[genename] == "":
			absdict[(strain)] = None
			continue

		elif genenames[genename] == "nan":
			absdict[(strain)] = None
			continue

		elif "refound" in genenames[genename]:
			absdict[(strain)] = None
			continue

		with open(f"{strain}.fasta", "w") as fasta:
			fasta.write(open(f"{gffpath}{strain}.gff", "r").read().split("##FASTA")[1])

		findex = Fasta(f"{strain}.fasta" , sequence_always_upper=True)

		with open(f"{gffpath}{strain}.gff", "r") as gff:
			if paralogs == True:
			
				for line in gff:

					if line.startswith("##FASTA"):
						break
						
					elif line.startswith("#"):
						continue

					splitline = line.rstrip("\n").split("\t")

					if splitline[2] == "CDS" and splitline[8].split(";")[0].split("=")[1] in names:
						strand = splitline[6]
						#negative strand + paralogs
						if strand == "-":
							gene = splitline[8].split(";")[0].split("=")[1]
							start = int(splitline[3])
							stop = int(splitline[4])
							genelength = stop-start
							contig = splitline[0]
							seq = str(findex[contig][start-1:stop].reverse.complement)
							offset = 0
							if shift:
								if not seq.startswith(shift[0]):
									offset = seq.find(shift[1])
									if offset == -1:
										print(f"could not find {shift[1]} in {gene} for {strain}")
										offset = 0

							if boundaries:
								if stop - boundaries[0] - offset > len(str(findex[contig])):
									downbound = len(str(findex[contig]))
									upbound = stop - boundaries[1]
									ups = stop - boundaries[0] - offset - len(str(findex[contig]))
								
								elif stop - boundaries[1] <= 0:
									downbound = stop - boundaries[0] - offset
									upbound = 1
									downs = -1 * (stop - boundaries[1])

								else:
									downbound = stop - boundaries[0] - offset
									upbound = stop - boundaries[1] - offset

							else:
								if stop - upstream - offset > len(str(findex[contig])):
									downbound = len(str(findex[contig]))
									upbound = start - downstream
									ups = stop - upstream - offset - len(str(findex[contig]))
								
								elif start - downstream <= 0:
									downbound = stop - upstream - offset
									upbound = 1
									downs = -1 * (start - downstream)
								
								else:
									downbound = stop - upstream - offset
									upbound = start - downstream
						
							seq = list(str(findex[contig][upbound-1:downbound].reverse.complement))
						
						else:
							#positive strand + paralogs
							gene = splitline[8].split(";")[0].split("=")[1]
							start = int(splitline[3])-1
							stop = int(splitline[4])-1
							strand = splitline[6]
							genelength = stop-start
							
							contig = splitline[0]
							seq = str(findex[contig][start-1:stop])
							offset = 0
							
							if shift:
								if not seq.startswith(shift[0]):
									offset = seq.find(shift[1])
									if offset == -1:
										print(f"could not find {shift[1]} in {gene} for {strain}")
										offset = 0
										
							if boundaries:
								if start + boundaries[0] + offset <= 0:
									upbound = 1
									downbound = start + boundaries[1]
									ups = -1 * (start + boundaries[0] + offset)
								
								elif start + boundaries[1] > len(str(findex[contig])):
									downbound = start + boundaries[0] + offset
									upbound = len(str(findex[contig]))
									downs = start + boundaries[1] - len(str(findex[contig]))

								else:
									upbound = start + boundaries[0] + offset
									downbound = start + boundaries[1] + offset
									
							else:
								if start + upstream + offset <= 0:
									downbound = stop + downstream
									upbound = 1 + offset
									ups = -1 * (start + upstream)
								
								elif stop + downstream > len(str(findex[contig])):
									downbound = len(str(findex[contig]))
									upbound = start + upstream + offset
									downs = stop + downstream - len(str(findex[contig]))
								
								else:
									downbound = stop + downstream
									upbound = start + upstream + offset

							seq = list(str(findex[contig][upbound-1:downbound]))
						
						if boundaries:
							if len(seq) < len(range(boundaries[0], boundaries[1])) and ups is not None:
								nanarr = np.empty(len(range(boundaries[0], boundaries[1])) - len(seq))
								nanarr[:] = np.nan
								seq = list(nanarr) + seq

							elif len(seq) < len(range(boundaries[0], boundaries[1])) and downs is not None:
								nanarr = np.empty(len(range(boundaries[0], boundaries[1])) - len(seq))
								nanarr[:] = np.nan
								seq = seq + list(nanarr)
						
						else:
							if ups is not None:
								nanarr = np.empty(ups)
								nanarr[:] = np.nan
								seq = list(nanarr) + seq

							elif downs is not None:
								nanarr = np.empty(downs)
								nanarr[:] = np.nan
								seq = seq + list(nanarr)
						
						if longestgene < len(seq):
							longestgene = len(seq)
							
						grapharray = []

						###converting nucleotide list into scalar array if nucleotides are not being displayed
						if not hybrid and not nucleotides:
							for letter in seq:
								if letter in base2int:
									grapharray.append(1)
								else:
									grapharray.append(0)

						###converting nucleotide list into float array as indices for nucleotides 
						else:
							for letter in seq:
								if letter in base2int:
									grapharray.append(base2int[letter])
								else:
									grapharray.append(float("nan"))

						seqdict[(strain, gene)] = grapharray
						shiftdict[(strain, gene)] = offset
						ups = None
						downs = None

			else:
				for line in gff:
					if line.startswith("##FASTA"):
						break
					elif line.startswith("#"):
						continue

					splitline = line.rstrip("\n").split("\t")

					if splitline[2] == "CDS" and splitline[8].split(";")[0].split("=")[1] == genenames[genename]:
						strand = splitline[6]
						if strand == "-":
							#negative strand without paralogs
							gene = splitline[8].split(";")[0].split("=")[1]
							start = int(splitline[3])
							stop = int(splitline[4])
							genelength = stop-start
							
							contig = splitline[0]
							seq = str(findex[contig][start-1:stop].reverse.complement)
							offset = 0
							if shift:
								if not seq.startswith(shift[0]):
									offset = seq.find(shift[1])
									if offset == -1:
										print(f"could not find {shift[1]} in {gene} for {strain}")
										offset = 0

							if boundaries:
								if stop - boundaries[0] - offset > len(str(findex[contig])):
									downbound = len(str(findex[contig]))
									upbound = stop - boundaries[1]
									ups = stop - boundaries[0] - offset - len(str(findex[contig]))
								
								elif stop - boundaries[1] <= 0:
									downbound = stop - boundaries[0] - offset
									upbound = 1
									downs = -1 * (stop - boundaries[1])

								else:
									downbound = stop - boundaries[0] - offset
									upbound = stop - boundaries[1] - offset

							else:
								if stop - upstream - offset > len(str(findex[contig])):
									downbound = len(str(findex[contig]))
									upbound = start - downstream
									ups = stop - upstream - offset - len(str(findex[contig]))
								
								elif start - downstream <= 0:
									downbound = stop - upstream - offset
									upbound = 1
									downs = -1 * (start - downstream)
								
								else:
									downbound = stop - upstream - offset
									upbound = start - downstream
						
							seq = list(str(findex[contig][upbound-1:downbound].reverse.complement))
								
						else:
							#positive strand without paralogs
							gene = splitline[8].split(";")[0].split("=")[1]
							start = int(splitline[3])-1
							stop = int(splitline[4])-1
							genelength = stop-start
							
							contig = splitline[0]
							seq = str(findex[contig][start-1:stop])

							offset = 0
							
							if shift:
								if not seq.startswith(shift[0]):
									offset = seq.find(shift[1])
									if offset == -1:
										print(f"could not find {shift[1]} in {gene} for {strain}")
										offset = 0
										
							if boundaries:
								if start + boundaries[0] + offset <= 0:
									upbound = 1
									downbound = start + boundaries[1]
									ups = -1 * (start + boundaries[0] + offset)
								
								elif start + boundaries[1] > len(str(findex[contig])):
									downbound = start + boundaries[0] + offset
									upbound = len(str(findex[contig]))
									downs = start + boundaries[1] - len(str(findex[contig]))

								else:
									upbound = start + boundaries[0] + offset
									downbound = start + boundaries[1] + offset
									
							else:
								if start + upstream + offset <= 0:
									downbound = stop + downstream
									upbound = 1 + offset
									ups = -1 * (start + upstream)
								
								elif stop + downstream > len(str(findex[contig])):
									downbound = len(str(findex[contig]))
									upbound = start + upstream + offset
									downs = stop + downstream - len(str(findex[contig]))
								
								else:
									downbound = stop + downstream
									upbound = start + upstream + offset

							seq = list(str(findex[contig][upbound-1:downbound]))
							#print(seq)
							
						if boundaries:
							if len(seq) < len(range(boundaries[0],boundaries[1])) and ups is not None:
								nanarr = np.empty(len(range(boundaries[0],boundaries[1])) - len(seq))
								nanarr[:] = np.nan
								seq = list(nanarr) + seq

							elif len(seq) < len(range(boundaries[0],boundaries[1])) and downs is not None:
								nanarr = np.empty(len(range(boundaries[0],boundaries[1])) - len(seq))
								nanarr[:] = np.nan
								seq = seq + list(nanarr)

						else:
							if ups is not None:
								nanarr = np.empty(ups)
								nanarr[:] = np.nan
								seq = list(nanarr) + seq

							elif downs is not None:
								nanarr = np.empty(downs)
								nanarr[:] = np.nan
								seq = seq + list(nanarr)
						
						if longestgene < len(seq):
							longestgene = len(seq)
						
						grapharray = []

						###converting nucleotide list into scalar array if nucleotides are not being displayed
						if not hybrid and not nucleotides:
							for letter in seq:
								if letter in base2int:
									grapharray.append(1)
								else:
									grapharray.append(0)

						###converting nucleotide list into float array as indices for nucleotides 
						else:
							for letter in seq:
								if letter in base2int:
									grapharray.append(base2int[letter])
								else:
									grapharray.append(float("nan"))

						seqdict[(strain, genenames[genename])] = grapharray
						shiftdict[(strain, genenames[genename])] = offset
						ups = None
						downs = None

		paralogs = False
	
	###here we append nans or zeros to genes, depending on whether nucleotides are being displayed. 
	###This step is needed to reach an equal length (necessary for plot creation).
	if boundaries:
		for geneinstance in seqdict:
			seqdict[geneinstance] = np.array(seqdict[geneinstance], dtype=float)
			
			if len(seqdict[geneinstance]) < len(range(boundaries[0],boundaries[0]+longestgene)):

				if not hybrid and not nucleotides:
					zeroarr = np.zeros(len(range(boundaries[0],boundaries[0]+longestgene)) - len(seqdict[geneinstance]))
					seqdict[geneinstance] = np.append(seqdict[geneinstance], zeroarr)

				else:
					nanarr = np.empty(len(range(boundaries[0],boundaries[0]+longestgene)) - len(seqdict[geneinstance]))
					nanarr[:] = np.nan
					seqdict[geneinstance] = np.append(seqdict[geneinstance], nanarr)

	else:
		for geneinstance in seqdict:
			seqdict[geneinstance] = np.array(seqdict[geneinstance], dtype=float)

			if len(seqdict[geneinstance]) < longestgene:

				if not hybrid and not nucleotides:
					zeroarr = np.zeros(longestgene - len(seqdict[geneinstance]))
					seqdict[geneinstance] = np.append(seqdict[geneinstance], zeroarr)

				else:
					nanarr = np.empty(longestgene - len(seqdict[geneinstance]))
					nanarr[:] = np.nan
					seqdict[geneinstance] = np.append(seqdict[geneinstance], nanarr)

	###here we start to gather the p-values from the sig_kmers_data.tsv
	pvaldict = {}
	if boundaries:
		for geneinstance in seqdict:
			pvaldict[geneinstance] = {}
			for x in range(boundaries[0], boundaries[0] + longestgene):
				pvaldict[geneinstance][x] = []
					
	else:
		for geneinstance in seqdict:
			pvaldict[geneinstance] = {}
			for x in range(upstream, longestgene + upstream):
				pvaldict[geneinstance][x] = []
		
	with open(sigkmers, "r") as sikm:
		sikm.readline()
		for line in sikm:
			splitline = line.rstrip().split("\t")
			if boundaries:
				if splitline[2] == clusterid and splitline[4] in genenames:
					#assigning vars for ease of understanding
					sigstrain = splitline[3]
					siggene = splitline[4]
					specshift = shiftdict[(sigstrain, siggene)]
					b0 = boundaries[0] + specshift
					b1 = boundaries[1] + specshift
					start = splitline[9]
					stop = splitline[10]
					
					if start == "" or stop == "":

						for x in range(boundaries[0],boundaries[1]):
							pvaldict[(sigstrain, siggene)][x].append(float(splitline[14]))
					
					elif b0 < int(start) and int(stop) < b1:#k-mer is entirely within boundaries

						for x in range(int(start)-specshift, int(stop)-specshift):
							pvaldict[(sigstrain, siggene)][x].append(float(splitline[14]))
					
					elif b0 < int(start) < b1 and b1 < int(stop):#k-mer starts in boundaries and ends outside

						for x in range(int(start)-specshift, boundaries[1]):
							pvaldict[(sigstrain, siggene)][x].append(float(splitline[14]))
					
					elif int(start) < b0 and b0 < int(stop) < b1:#k-mer starts outside of boundaries and ends inside

						for x in range(boundaries[0], int(stop)-specshift):
							pvaldict[(sigstrain, siggene)][x].append(float(splitline[14]))
							
					elif int(start) < b0 and b1 < int(stop):#boundaries lie within longer k-mer
						
						for x in range(boundaries[0],boundaries[1]):
							pvaldict[(sigstrain, siggene)][x].append(float(splitline[14]))
					
					else:

						continue

			else:
				if splitline[2] == clusterid and splitline[4] in genenames:
					sigstrain = splitline[3]
					siggene = splitline[4]
					#shift as identified from the nucleotide sequence for a specific (strain, genename)-tuple
					specshift = shiftdict[(sigstrain, siggene)]
					start = splitline[9]
					stop = splitline[10]
					contigstart = int(splitline[7])
					contigstop = int(splitline[8])

					#triggers if the entire gene has been flagged as significantly associated
					if start == "" or stop == "":
						for x in range(upstream,longestgene+upstream-specshift):
							pvaldict[(sigstrain, siggene)][x].append(float(splitline[14]))

					#triggers for significantly associated k-mers
					else:
						for x in range(int(start)-specshift,int(stop)-specshift):
							pvaldict[(sigstrain, siggene)][x].append(float(splitline[14]))

	###calculate threshold
	with open(htpfile, "r") as htp:
		linecount = 0
		next(htp)
		for line in htp:
			linecount += 1
			
		#remove count for trailing "\n"
		linecount -= 1

	threshold = -np.log10(alpha / float(linecount))

	###using partial to make the following calculation of scalar values more readable (or less unreadable...)
	scal_func = partial(calc_scalar, 
							maxval = maxscalar, 
							thresh = threshold)

	###calculate mean p-value for each position in the gene and subsequently calculating the -log10 as well as the scalar value (used to display the opacity)
	if hybrid:
		for geneinstance in pvaldict:
			pvaldict[geneinstance] = {x: (scal_func(-np.log10(mean(pvaldict[geneinstance][x])))) if len(pvaldict[geneinstance][x]) > 0 else 0.1 for x in pvaldict[geneinstance]}
			pvaldict[geneinstance] = np.array([*zip(*pvaldict[geneinstance].items())][1], dtype = float)
	
	elif nucleotides:
		pass

	else:###calculate mean p-value for each position in the gene and subsequently calculating the -log10
		for geneinstance in pvaldict:
			pvaldict[geneinstance] = {x: (-np.log10(mean(pvaldict[geneinstance][x]))) if len(pvaldict[geneinstance][x]) > 0 else "NaN" for x in pvaldict[geneinstance]}
			pvaldict[geneinstance] = np.array([*zip(*pvaldict[geneinstance].items())][1], dtype = float)

	###creating arrays for strains that are not present within the gene cluster and setting the order in which strains are displayed in the plot
	nanarr = np.empty(longestgene)
	zeroarray = np.zeros(longestgene)
	zeroabsdict = {key: zeroarray for key in absdict.keys()}
	nanarr[:] = np.nan
	absdict = {key: nanarr for key in absdict.keys()}

	###create the x- and y-values for the plot/ indices and column names for the data frame
	firstlist = []
	lastlist = []
			
	if phenotype:
		for geneinstance in seqdict.keys():
			if geneinstance[0] in positivelist:
				firstlist.append(geneinstance)
			else:
				lastlist.append(geneinstance)
		yvals = firstlist + lastlist + list(absdict.keys())

	else:
		priodict = {}
		lastlist = []
		for geneinstance in seqdict:
			if geneinstance[0] in orderdict:
				if orderdict[geneinstance[0]] in priodict:
					priodict[orderdict[geneinstance[0]]].append(geneinstance)
				else:
					priodict[orderdict[geneinstance[0]]] = [geneinstance,]
			else:
				lastlist.append(geneinstance)

		lastlist = lastlist + list(absdict.keys())

		firstlist

		for phenotypevalue in sorted(priodict):
			firstlist.extend(priodict[phenotypevalue])

		yvals = firstlist + lastlist

	if boundaries:
		xvals = [*range(boundaries[0], boundaries[0] + longestgene)]

	else:
		xvals = [*range(upstream, longestgene + upstream)]

	###transfer pvaldict into pandas dataframe
	pvaldf = pd.DataFrame(index=yvals, columns=xvals)

	for geneinstance in pvaldict:
		pvaldf.loc[[geneinstance],] = pvaldict[geneinstance]

	for strain in zeroabsdict:
		pvaldf.loc[strain,] = zeroabsdict[strain]

	pvalarray = pvaldf.to_numpy(dtype = "float")

	###transfer seqdict into pandas dataframe
	seqdf = pd.DataFrame(index=yvals, columns=xvals)

	for geneinstance in seqdict:
		seqdf.loc[[geneinstance],] = seqdict[geneinstance]

	if not hybrid and not nucleotides: 
		for strain in zeroabsdict:
			seqdf.loc[strain,] = zeroabsdict[strain]

	else:
		for strain in absdict:
			seqdf.loc[strain,] = absdict[strain]
	
	seqarray = seqdf.to_numpy(dtype = "float")

	###create color map
	if hybrid or nucleotides:
		base_colors = list(sns.color_palette('tab20', 4))
		cmap = colors.LinearSegmentedColormap.from_list('nucleotides', base_colors, 4)
		cmap.set_bad('xkcd:grey')

	else:
		cmap = copy.copy(plt.get_cmap('viridis'))
		cmap.set_bad('xkcd:grey')
		cmap.set_under('xkcd:light grey')
	
	fig, ax = plt.subplots()

	###create the heatmap
	if hybrid:
		im = ax.imshow(seqarray, cmap = cmap, aspect="auto", alpha=pvalarray)
	elif nucleotides:
		im = ax.imshow(seqarray, cmap = cmap, aspect="auto")
	else:
		im = ax.imshow(pvalarray, cmap = cmap, vmin = threshold, vmax = maxscalar, aspect="auto", alpha=seqarray)

	###modify ticks and tick labels
	if (hybrid or nucleotides) and not display_nucs: 
		values = np.unique(seqarray.ravel())
		colors = [ im.cmap(im.norm(value)) for value in values]
		patches = [ mpatches.Patch(color=colors[i], label=f"{float2legend[str(values[i])]}".format(l=values[i]) ) for i in range(len(values)) ]
		#put those patched as legend-handles into the legend
		plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
	elif (hybrid or nucleotides) and display_nucs:
		for (j,i),label in np.ndenumerate(seqarray):

			ax.text(i,j,float2legend[str(label)],ha='center',va='center')

	###set tick labels
	ax.set_xticks(np.arange(len(xvals)), labels=xvals)

	if boundaries is None:
		n = 100	 # Keeps every 100th label
		[l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % n != 0]

		ax.axes.get_yaxis().set_visible(False)
	else:
		n = 20	 # Keeps every 20th label
		[l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % n != 0]
	
		ax.axes.get_yaxis().set_visible(False)

	ax.tick_params(width = 1, labelsize="x-small")
	mpl.pyplot.xticks(rotation=90)

	###draw lines to show upstream/downstream region
	if not boundaries:
		xcoords = [-upstream, longestgene+upstream]
		for xc in xcoords:
			ax.axvline(x=xc, lw = 0.5, color = "black")

	###draw line to show phenotype positive strains
	if strainmarker:
		ax.axhline(y=len(positivelist)-0.5, lw = 0.5, color = "black")
	
	###colorbar for p-values only
	if not hybrid and not nucleotides:
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		coba = plt.colorbar(im, cax=cax)
		coba.set_label("-log10 p-value")

	###set figure size
	mpl.pyplot.rcParams["figure.figsize"] = [10, 9]

	###labeling of axes and title
	ax.set_title(f"significant k-mers {clusterid}")
	ax.set_xlabel("position")
	ax.set_ylabel("strains")

	mpl.pyplot.style.use("ggplot")

	###get current figure. is this necessary? try deleting it if I find time
	fig1 = plt.gcf()
	
	###name file according to options used
	if hybrid:
		savename = "hybrid_heatmap"
	elif nucleotides:
		savename = "nucleotides_plot"
	else:
		savename = "p-value_plot"

	if shift:
		savename = savename + "_shifted"

	###save file
	if boundaries:
		plt.savefig(f"{outpath}{clusterid}_{savename}_{boundaries[0]}-{boundaries[1]}.png", dpi=1400, bbox_inches = "tight")
	else:
		plt.savefig(f"{outpath}{clusterid}_{savename}.png", dpi=1400, bbox_inches = "tight")

	###is this necessary? try deleting it if I find time
	plt.figure().clear()
	plt.close()
	plt.cla()
	plt.clf()

	###remove temporary fasta and fai files in directory
	for file in os.listdir(outpath):
		if file.endswith(".fasta") or file.endswith(".fai") or file.startswith("temp_"):
			os.remove(f"{outpath}{file}")
