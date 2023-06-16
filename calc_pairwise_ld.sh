#!/bin/bash
#
# Usage:
#   calc_pairwise_ld.sh <indir> <outdir> <ncpus_ngsld>
#
# Description:
#   Calculate pairwise linkage disequilibrium (LD) with ngsLD.
#
# Requirements:
#   ngsld

# path to ngsLD
ngsld=${HOME}/software/ngsLD/ngsLD

# number of samples
n_ind=$(cat $1/bamlist | wc -l)
# find MAF file
maf=$(find $1 -name '*.mafs.gz')
# create file with SNP positions
zcat ${maf} | cut -f 1,2 | tail -n +2 > $2/snps.pos
# count the number of SNPs
n_sites=$(cat $2/snps.pos | wc -l)
# find genotype likelihoods (GL) file
gl=$(find $1 -name '*.beagle.gz')

# calculate LD for all SNP pairs up to 500 kbb apart
${ngsld} \
  --geno ${gl} \
  --probs \
  --n_ind ${n_ind} \
  --n_sites ${n_sites} \
  --pos $2/snps.pos \
  --max_kb_dist 500 \
  --out $2/snps.ld \
  --n_threads $3 \
  &> $2/ngsld.log

# sample 0.05% of all estimated SNPs pairs for fitting LD decay curves later on
awk 'rand()<0.0005' $2/snps.ld > $2/snps.ld.sampled