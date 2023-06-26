#!/bin/bash
argumentarray=()
##############single argument changes
#basic
panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --output comp_data/panfeed_out/basic &> comp_data/logs/basic.log
argumentarray+=("basic")
#no k-mers logged
panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --output comp_data/panfeed_out/nolog &> comp_data/logs/nolog.log
argumentarray+=("nolog")
#upstream
panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --upstream 100 --downstream 0 --output comp_data/panfeed_out/upstream &> comp_data/logs/upstream.log
argumentarray+=("upstream")
#downstream
panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --upstream 0 --downstream 100 --output comp_data/panfeed_out/downstream &> comp_data/logs/downstream.log
argumentarray+=("downstream")
#up-/downstream
panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --upstream 100 --downstream 100 --output comp_data/panfeed_out/updownstream &> comp_data/logs/updownstream.log
argumentarray+=("updownstream")
#up-/downstream
panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --upstream 100 --downstream 100 --downstream-start-codon --output comp_data/panfeed_out/downstart &> comp_data/logs/downstart.log
argumentarray+=("downstart")
#non-canonical
panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --non-canonical --output comp_data/panfeed_out/noncanonical &> comp_data/logs/noncanonical.log
argumentarray+=("noncanonical")
#no filter
panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --no-filter --output comp_data/panfeed_out/nofilter &> comp_data/logs/nofilter.log
argumentarray+=("nofilter")
#high maf
panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --maf 0.1 --output comp_data/panfeed_out/highmaf &> comp_data/logs/highmaf.log
argumentarray+=("highmaf")
#consider missing
panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --consider-missing --output comp_data/panfeed_out/considermissing &> comp_data/logs/considermissing.log
argumentarray+=("considermissing")
#compressed
#panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --compress --output comp_data/panfeed_out/compress &> comp_data/logs/compress.log
#argumentarray+="compress"
#multiple files
#panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --multiple-files --output comp_data/panfeed_out/multiplefiles &> comp_data/logs/multiplefiles.log
#argumentarray+="multiplefiles"
#multiple cores
#panfeed --gff ./gffs/ --presence-absence ./gene_presence_absence.csv --genes target_clusters.txt --targets stroiall.txt --cores 3 --output comp_data/panfeed_out/multiplecores &> comp_data/logs/multiplecores.log
#argumentarray+="multiplecores"
##############multiple argument changes
##############compare output files
for argument in "${argumentarray[@]}";
do
  echo "Comparing results and logs to baseline files "$argument;
  awk -F "-" '{$1=""; print $0}' comp_data/logs/$argument.log > tmp_log.log && mv tmp_log.log comp_data/logs/$argument.log
  if cmp -s comp_data/logs/$argument.log comp_data/logs/baseline_$argument.log; then
  	echo $argument" logs are the same";
  else
  	echo $argument" logs are different";
  fi

  if cmp -s comp_data/panfeed_out/$argument/kmers.tsv comp_data/panfeed_out/baseline_$argument/kmers.tsv \
  	&& cmp -s comp_data/panfeed_out/$argument/kmers_to_hashes.tsv comp_data/panfeed_out/baseline_$argument/kmers_to_hashes.tsv \
	&& cmp -s comp_data/panfeed_out/$argument/hashes_to_patterns.tsv comp_data/panfeed_out/baseline_$argument/hashes_to_patterns.tsv; then
	echo $argument" files are the same";
  else
  	echo $argument" files are different";
  fi
done
