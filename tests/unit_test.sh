#!/bin/bash
argumentarray=()

green="\033[1;32m"
red="\033[1;31m"
reset="\033[0m"

die () {
        echo -e $red"############"$reset
	echo -e $red$1$reset
	echo -e $red"Test failed!"$reset
	echo -e $red"############"$reset
	exit 1
}
##############single argument changes
#basic
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --output test_files/comp_data/panfeed_out/basic &> test_files/comp_data/logs/basic.log || die "Basic run failed"
argumentarray+=("basic")
#cores
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --cores 4 --output test_files/comp_data/panfeed_out/cores &> test_files/comp_data/logs/cores.log || die "Multiple cores run failed"
argumentarray+=("cores")
#no k-mers logged
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --output test_files/comp_data/panfeed_out/nolog &> test_files/comp_data/logs/nolog.log || die "No k-mers logged run failed"
argumentarray+=("nolog")
#upstream
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --upstream 100 --downstream 0 --output test_files/comp_data/panfeed_out/upstream &> test_files/comp_data/logs/upstream.log || die "Upstream run failed"
argumentarray+=("upstream")
#downstream
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --upstream 0 --downstream 100 --output test_files/comp_data/panfeed_out/downstream &> test_files/comp_data/logs/downstream.log || die "Downstream run failed"
argumentarray+=("downstream")
#up-/downstream
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --upstream 100 --downstream 100 --output test_files/comp_data/panfeed_out/updownstream &> test_files/comp_data/logs/updownstream.log || die "Up-/downstream run failed"
argumentarray+=("updownstream")
#up-/downstream
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --upstream 100 --downstream 100 --downstream-start-codon --output test_files/comp_data/panfeed_out/downstart &> test_files/comp_data/logs/downstart.log || die "Downstream start codon run failed"
argumentarray+=("downstart")
#non-canonical
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --non-canonical --output test_files/comp_data/panfeed_out/noncanonical &> test_files/comp_data/logs/noncanonical.log || die "Non-canonical run failed"
argumentarray+=("noncanonical") 
#no filter
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --no-filter --output test_files/comp_data/panfeed_out/nofilter &> test_files/comp_data/logs/nofilter.log || die "No filter run failed"
argumentarray+=("nofilter")
#high maf
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --maf 0.1 --output test_files/comp_data/panfeed_out/highmaf &> test_files/comp_data/logs/highmaf.log || die "High MAF run failed"
argumentarray+=("highmaf")
#consider missing
python ../panfeed-runner.py --gff test_files/gffs/ --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --consider-missing --output test_files/comp_data/panfeed_out/considermissing &> test_files/comp_data/logs/considermissing.log || die "Consider missing run failed"
argumentarray+=("considermissing")
#fileoffiles
python ../panfeed-runner.py --gff test_files/input_gffs.txt --fasta test_files/input_fastas.txt --presence-absence test_files/gene_presence_absence.csv --targets test_files/stroi.txt --output test_files/comp_data/panfeed_out/fileoffiles &> test_files/comp_data/logs/fileoffiles.log || die "File of files run failed"
argumentarray+=("fileoffiles")
##############multiple argument changes
##############compare output files
for argument in "${argumentarray[@]}";
do
  echo "Comparing results and logs to baseline files "$argument;
  awk -F "-" '{$1=""; print $0}' test_files/comp_data/logs/$argument.log > tmp_log.log && mv tmp_log.log test_files/comp_data/logs/$argument.log
  if cmp -s test_files/comp_data/logs/$argument.log test_files/comp_data/logs/baseline_$argument.log; then
  	echo $argument" logs are the same";
  else
  	echo $argument" logs are different";
  fi
  if cmp -s test_files/comp_data/panfeed_out/$argument/kmers.tsv test_files/comp_data/panfeed_out/baseline_$argument/kmers.tsv \
  	&& cmp -s test_files/comp_data/panfeed_out/$argument/kmers_to_hashes.tsv test_files/comp_data/panfeed_out/baseline_$argument/kmers_to_hashes.tsv \
	&& cmp -s test_files/comp_data/panfeed_out/$argument/hashes_to_patterns.tsv test_files/comp_data/panfeed_out/baseline_$argument/hashes_to_patterns.tsv; then
	echo $argument" files are the same";
  else
  	echo $argument" files are different";
  fi
done
