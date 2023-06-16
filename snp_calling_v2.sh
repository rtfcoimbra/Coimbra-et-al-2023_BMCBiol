#!/bin/bash
#
# Usage:
#   snp_calling_v2.sh <cluster|all> <ref.fa> <dir_bams> <site_depth.stats> <dir_pruned_sites> <outdir> <njobs>
#
# Description:
#   Call and filter SNPs in parallel by chromosome with ANGSD.
#
# Requirements:
#   angsd
#   parallel
#   python3
#
# Important:
#   Each job uses up to 4 CPU threads.

mode=$1
ref=$2
dir_bams=$3
depth_stats=$4
dir_pruned_sites=$5
outdir=$6
njobs=$7

# check running mode
if [[ ${mode} == 'cluster' ]]; then
  # create list of input BAMs for a cluster
  ls -1 -v ${dir_bams}/*.clean.bam > ${outdir}/bamlist
elif [[ ${mode} == 'all' ]]; then
  # create list of input BAMs for all individuals
  ls -1 -v ${dir_bams}/{masai,nubian,reticulated}/*.clean.bam > ${outdir}/bamlist
else
  echo "Unrecognized argument: ${mode} ... expected 'cluster' or 'all'"
  exit 1
fi

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
filter_snps="-sb_pval 1e-6 -hwe_pval 1e-6 -hetbias_pval 1e-6 -SNP_pval 1e-6"
# set angsd tasks
todo="-doSnpStat 1 -doHWE 1 -doCounts 1 -doMajorMinor 1 -doMaf 2 -doPost 1 -doGeno 8 -doGlf 2"

# check running mode
if [[ ${mode} == 'cluster' ]]; then
  # SNP calling per cluster with min. MAF and SNP filters
  parallel -j ${njobs} \
    "angsd -b ${outdir}/bamlist -ref ${ref} -r {} ${filter_reads} ${filter_sites} -minMaf 0.05 ${filter_snps} ${todo} -GL 1 -P 4 -out ${outdir}/{}.snps &> ${outdir}/{}.angsd.log" \
    ::: chr{1..14}
elif [[ ${mode} == 'all' ]]; then
  # SNP calling for all individuals without min. MAF and SNP filters
  parallel -j ${njobs} \
    "angsd -b ${outdir}/bamlist -ref ${ref} -r {} -sites ${dir_pruned_sites}/combined.snps.ld.pruned.{} ${filter_reads} ${filter_sites} ${todo} -GL 1 -P 4 -out ${outdir}/{}.snps &> ${outdir}/{}.angsd.log" \
    ::: chr{1..14}
else
  echo "Unrecognized argument: ${mode} ... expected 'cluster' or 'all'"
  exit 1
fi

# concatenate output files
gunzip ${outdir}/*.gz
for type in beagle mafs hwe snpStat; do
  head -1 ${outdir}/chr1.snps.${type} > ${outdir}/snps.${type}
done
cat ${outdir}/chr{1..14}.snps.beagle | sed '/marker/d' >> ${outdir}/snps.beagle
cat ${outdir}/chr{1..14}.snps.mafs | sed '/chromo/d' >> ${outdir}/snps.mafs
cat ${outdir}/chr{1..14}.snps.hwe | sed '/Chromo/d' >> ${outdir}/snps.hwe
cat ${outdir}/chr{1..14}.snps.snpStat | sed '/Chromo/d' >> ${outdir}/snps.snpStat
cat ${outdir}/chr{1..14}.snps.geno >> ${outdir}/snps.geno
for type in beagle geno mafs hwe snpStat; do
  gzip ${outdir}/snps.${type}
done

# generate optional sites file for bcftools
#zcat ${outdir}/snps.beagle.gz \
#  | awk 'NR>1 { print $1 }' \
#  | sed 's/_/\t/' > ${outdir}/snps.sites
#parallel -j ${njobs} \
#  "grep -P '{}\t' ${outdir}/snps.sites > ${outdir}/snps.sites.{}" \
#  ::: chr{1..14}

# clean directory
rm ${outdir}/chr{1..14}.snps.{beagle,geno,mafs,hwe,snpStat} #${outdir}/snps.sites