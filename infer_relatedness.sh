#!/bin/bash
#
# Usage:
#   infer_relatedness.sh <path_to_maf_beagle> <indir> <outdir> <k> <min_maf> <ncpus>
#
# Description:
#   Estimate individual ancestries and allele frequencies for a given K value
#   with 100 replicates in NGSadmix, find the highest likelihood run, and infer
#   relatedness among individuals based on it with NGSremix.
#
# Requirements:
#   maf_beagle (https://github.com/KHanghoj/maf_beagle_ngsadmix)
#   ngsadmix
#   ngsremix

# set variables from positional arguments
maf_beagle=$1
indir=$2
outdir=$3
k=$4
min_maf=$5
ncpus=$6

# find genotype likelihoods file
gl=$(find ${indir} -name '*.beagle.gz')

# create ngsadmix output directory
outdir_ngsadmix=${outdir}/ngsadmix_k${k}
mkdir -p ${outdir_ngsadmix}

# set ngsadmix output file name and directory
fout=${outdir_ngsadmix}/$(basename ${gl%.beagle.gz})

# estimate individual ancestries and allele frequencies
# for a given K with 100 replicates
for run in {1..100}; do
  NGSadmix \
    -likes ${gl} \
    -K ${k} \
    -maxiter 5000 \
    -minMaf ${min_maf} \
    -o ${fout}.k${k}_r${run} \
    -P ${ncpus}
done

# list ngsadmix log files
logs=$(ls -v ${fout}.k${k}_r*.log)

# find the highest likelihood run
best_run=$(( for log in ${logs}; do grep -Po 'like=\K[^ ]+' ${log}; done ) | sort -gr | head -1 | grep -f - ${logs} | grep -Po '_r\K[0-9]+')

# extract MAF filtered sites in ngsadmix from genotype likelihoods file
# (https://github.com/KHanghoj/NGSremix/issues/3#issuecomment-1135490280)
zcat ${gl} | ${maf_beagle} ${min_maf} 2> ${outdir}/maf_beagle.log | gzip -c > ${outdir}/snps.minmaf.beagle.gz

# estimate relatedness
NGSremix \
  -beagle ${outdir}/snps.minmaf.beagle.gz \
  -qname ${fout}.k${k}_r${best_run}.qopt \
  -fname ${fout}.k${k}_r${best_run}.fopt.gz \
  -o ${outdir}/ngsremix \
  -P ${ncpus}