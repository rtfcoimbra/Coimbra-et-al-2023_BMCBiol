#!/bin/bash
#
# Usage:
#   process_pruned_snps.sh <ld_pruning_dir> <njobs>
#
# Description:
#   Combine pruned SNPs of all species while removing duplicated sites, split
#   them into different files by chromosome, and index the resulting files for
#   ANGSD.
#
# Requirements:
#   angsd
#   parallel

# set variables from positional arguments
ld_pruning_dir=$1
njobs=$2

# create array of input files
fins=(${ld_pruning_dir}/{masai,nubian,reticulated}/snps.ld.pruned)

# set output file path
fout=${ld_pruning_dir}/combined.snps.ld.pruned

# combine pruned SNP sites of all species and remove duplicates
cat ${fins[@]} | sort -V | uniq > ${fout}

# split sites into different files by chromosome
parallel -j ${njobs} "grep -P {}_ ${fout} > ${fout}.{}" ::: chr{1..14}
sed -i 's/_/\t/' ${fout}.chr*

# index sites for angsd
for f in ${fout}.chr*; do angsd sites index ${f}; done