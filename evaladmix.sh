#!/bin/bash
#
# Usage:
#   evaladmix.sh <gl_dir> <admix_dir> <outdir> <k_min> <k_max> <min_maf> <ncpus>
#
# Description:
#   Calculate a matrix of pairwise correlation of residuals between
#   individuals using evalAdmix to test the goodness of fit of the data
#   to an admixture model.
#
# Requirements:
#   evaladmix

# find genotype likelihoods file
gl=$(find $1 -name '*.beagle.gz')

# set variables from positional arguments
admix_dir=$2
outdir=$3
k_min=$4
k_max=$5
min_maf=$6
ncpus=$7

# iterate over a range of K values
for k in $(seq ${k_min} ${k_max}); do
  # list all log files
  k_logs=$(ls -v ${admix_dir}/*.k${k}_r*.log)
  # find the highest likelihood run
  best_run=$(head -1 ${admix_dir}/likelihoods_k${k} | grep -f - ${k_logs} | grep -Po '_r\K[0-9]+')
  # get ngsadmix files
  fopt=${admix_dir}/snps.k${k}_r${best_run}.fopt.gz
  qopt=${admix_dir}/snps.k${k}_r${best_run}.qopt
  # output file
  fout=${outdir}/k${k}_r${best_run}.corres.txt
  # calculate a pairwise correlation of residuals matrix between individuals
  evalAdmix -beagle ${gl} -fname ${fopt} -qname ${qopt} -minMaf ${min_maf} -o ${fout} -P ${ncpus}
done