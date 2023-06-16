#!/bin/bash
#
# Usage:
#   estimate_fbranch.sh </path/dsuite> <vcf> <tree.nwk> <sample_sets.txt> <outdir>
#
# Description:
#   Calculate Patterson’s D, f4-ratio, and f-branch statistics for all
#   possible trios of population/species with Dsuite.
#
# Requirements:
#   dsuite (including 'dtools.py')
#   python3

dsuite=$1
vcf=$2
nwk_tree=$3
sample_sets=$4
outdir=$5

out_prefix=${outdir}/dsuite
run_name=$(basename ${nwk_tree} .rooted.nwk)
fout=${out_prefix}_${run_name}

# calculate the D and f4-ratio statistics for all possible trios of populations/species
Dsuite Dtrios -c -g -k 100 -o ${out_prefix} -n ${run_name} -t ${nwk_tree} ${vcf} ${sample_sets}

# disentangle correlated f4-ratio estimates and assign evidence for
# introgression to specific, possibly internal, branches on a phylogeny
# given that they can be tested under a ((P1, P2) P3, Outgroup) topology
Dsuite Fbranch ${nwk_tree} ${fout}_tree.txt > ${fout}_Fbranch.txt

# plot f-branch statistic as .png and .svg
python3 ${dsuite}/utils/dtools.py --dpi 300 ${fout}_Fbranch.txt ${nwk_tree}