#!/bin/bash
#
# Usage:
#   ld_pruning.sh <indir_gl> <indir_ld> <max_kb_dist> <min_weight> <njobs>
#
# Description:
#   Prune linked sites with 'prune_graph.pl' and extract unlinked sites from
#   ANGSD's genotype likelihoods file.
#
# Requirements:
#   parallel
#   perl
#   ngsld (incl. prune_graph.pl)

# path to 'prune_graph.pl'
prune_graph=${HOME}/software/ngsLD/scripts/prune_graph.pl
# path to local perl library
perl_lib=${HOME}/perl5/lib/perl5

# find genotype likelihoods (GL) file
gl=$(find $1 -name '*.beagle.gz')

# split ngsLD output into separate files by scaffold
parallel -j $5 "grep -P {}: $2/snps.ld > $2/snps.ld.{}" ::: chr{1..14}

# prune linked sites by scaffolds in parallel
parallel -j $5 \
  "perl -I ${perl_lib} ${prune_graph} --in_file {} --max_kb_dist $3 --min_weight $4 --out {}.pruned &> {}.prune_graph.log" \
  ::: $2/snps.ld.chr{1..14}

# concatenate and sort unlinked sites into a single file
cat $(ls $2/snps.ld.chr*.pruned) | sort -V | sed 's/:/_/' > $2/snps.ld.pruned
# concatenate 'prune_graph.pl' log files
cat $(ls -v $2/snps.ld.chr*.prune_graph.log) > $2/prune_graph.log
# remove intermediate files
rm $2/snps.ld.chr*

# set output file name
pruned_gl=$2/snps.ld_pruned.beagle
# extract SNPs matching the '*.ld.pruned' file from the GL file
zcat ${gl} | grep -F -f $2/snps.ld.pruned > ${pruned_gl}.tmp
# get linked SNPs inadvertently extracted due to partial string matching
awk '{ print $1 }' ${pruned_gl}.tmp \
  | diff $2/snps.ld.pruned - \
  | grep -Eo 'chr.+' > $2/sites2remove
# extract the header of the original GL file
zcat ${gl} | head -1 > ${pruned_gl}
# extract only the unlinked SNPs from the GL (.tmp) file
grep -v -F -f $2/sites2remove ${pruned_gl}.tmp >> ${pruned_gl}
# get unlinked SNPs inadvertently removed due to partial string matching
awk '{ print $1 }' ${pruned_gl} \
  | diff $2/snps.ld.pruned - \
  | grep -Eo 'chr.+' > $2/sites2recover
grep -F -f $2/sites2recover ${pruned_gl}.tmp >> ${pruned_gl}
# sort lines by SNP position respecting the header line and compress file
sed 's/marker/chr0/' ${pruned_gl} \
  | sort -V -k 1,1 \
  | sed 's/chr0/marker/' \
  | gzip > ${pruned_gl}.gz
# remove intermediate files
rm $2/sites2remove $2/sites2recover ${pruned_gl} ${pruned_gl}.tmp