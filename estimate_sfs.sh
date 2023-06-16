#!/bin/bash
#
# Usage:
#   estimate_sfs.sh <ref.fa> <anc.fa> <dir_bams> <site_depth.stats> <outdir> <ncpus>
#
# Description:
#   Estimate the unfolded site frequency spectrum (SFS) with ANGSD and realSFS.
#
# Requirements:
#   angsd
#   python3

ref=$1
anc=$2
dir_bams=$3
depth_stats=$4
outdir=$5
outname=$(basename ${dir_bams})
ncpus=$6

# create list of input BAMs
ls -1 -v ${dir_bams}/*.clean.bam > ${outdir}/bamlist

# count the number of individuals
n_ind=$(cat ${outdir}/bamlist | wc -l)
# set minimum number of individuals per site
min_ind=$(python -c "print(f'{round(${n_ind} * 0.9)}')")
# set minimum and maximum depth thresholds
min_dp=$(grep -Po '\(MEDIAN \- 5 \* MAD\): \K\-*\d+' ${depth_stats})
max_dp=$(grep -Po '\(MEDIAN \+ 5 \* MAD\): \K\d+' ${depth_stats})
# if minimum depth threshold is negative make it 1
if [[ ${min_dp} -lt 0 ]]; then min_dp=1; fi

# set read filters
filter_reads="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -baq 2"
# set site filters
filter_sites="-minMapQ 30 -minQ 30 -minInd ${min_ind} -setMinDepth ${min_dp} -setMaxDepth ${max_dp} -skipTriallelic 1"
# set SNP filters
filter_snps="-sb_pval 1e-6 -hetbias_pval 1e-6"
# set angsd tasks
todo="-doSnpStat 1 -doHWE 1 -doCounts 1 -doMajorMinor 1 -doMaf 2 -doPost 1 -doGeno 8 -doSaf 1"

# estimate site allele frequency (SAF) likelihoods per site
angsd \
  -b ${outdir}/bamlist \
  -ref ${ref} \
  -anc ${anc} \
  ${filter_reads} \
  ${filter_sites} \
  ${filter_snps} \
  ${todo} \
  -GL 1 \
  -P 4 \
  -out ${outdir}/${outname} \
  &> ${outdir}/angsd.log

# estimate the SFS
realSFS \
  ${outdir}/${outname}.saf.idx \
  -cores ${ncpus} \
  -maxiter 5000 \
  > ${outdir}/${outname}.sfs \
  2> ${outdir}/realSFS.log