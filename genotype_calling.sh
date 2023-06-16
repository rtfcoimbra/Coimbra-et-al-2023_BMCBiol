#!/bin/bash
#
# Usage:
#   genotype_calling.sh <phylo|mig> <ref.fa> <dir_bams> <sample_sets.txt> <site_depth.stats> <dir_sites> <outdir> <njobs>
#
# Description:
#   Call genotypes with bcftools.
#
# Requirements:
#   bcftools
#   parallel

mode=$1
ref=$2
dir_bams=$3
sample_sets=$4
depth_stats=$5
dir_sites=$6
outdir=$7
njobs=$8

# set shell extglob
shopt -s extglob

# check running mode
if [[ ${mode} == 'phylo' ]]; then
  # create list of giraffe and okapi input BAMs
  ls -1 -v ${dir_bams}/{GNP01,GNP04,GNP05,SNR2,ZNP01,ETH1,ETH2,ETH3,GF084,GF089,GF090,GF094,GF101,GF105,MF06,MF22,MF24,GF220,GF221,GF224,GF258,GF260,GF261,GF263,GF264,GF267,GF268,GF269,GF270,GF271,GF272,GF081,GF095,WA720,WA733,WA746,WA806,WA808,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110,BNP02,KKR01,KKR08,MTNP09,V23,SUN3,GF289,GF211,GF213,GF214,GF290,GF227,GF228,GF229,GF230,GF232,GF234,GF235,GF236,GF237,ISC04,ISC08,GF292,GF295,GF276,GF277,GF278,GF281,GF283,GF286,GF288,GF187,GF201,GF203,GF206,GF218,GF223,LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,GF155,GF080,GF048,GF054,GF070,GF074,GF076,GF250,GF253,GF068,GF146,GF242,GF245,GF248,GF249,GF049,GF077,GF246,MA1,GF005,GF007,GF011,GF016,GF023,GF034,GF035,GF036,GF238,GF239,GF044,GF065,GF144,SGR01,SGR05,SGR07,SGR13,SGR14,GF059,WOAK}.clean.bam
elif [[ ${mode} == 'mig' ]]; then
  # create list of giraffe input BAMs
  ls -1 -v ${dir_bams}/{nubian,reticulated,masai}/!(GF085|GF097|GF109|GF114|GF115|GF117|GF126|GF130|GF132|GF137|GF156|GF159|GF164|GF165|GF168|GF182|GF189|GF193|GF194|GF233|ISC01).clean.bam > ${outdir}/bamlist
else
  echo "Unrecognized argument: ${mode} ... expected 'phylo' or 'mig'"
  exit 1
fi

# unset shell extglob
shopt -u extglob

# set minimum and maximum depth thresholds
min_dp=$(grep -Po '\(MEDIAN \- 5 \* MAD\): \K\-*\d+' ${depth_stats})
max_dp=$(grep -Po '\(MEDIAN \+ 5 \* MAD\): \K\d+' ${depth_stats})
# if minimum depth threshold is negative make it 1
if [[ ${min_dp} -lt 0 ]]; then min_dp=1; fi


if [[ ${mode} == 'phylo' ]]; then
  # joint genotype calling
  parallel -j ${njobs} \
    "bcftools mpileup -b ${outdir}/bamlist -D -f ${ref} -q 30 -Q 30 -r {} --ns UNMAP,SECONDARY,QCFAIL,DUP --lu PROPER_PAIR -a AD -Ou \
    | bcftools call -mv -G ${sample_sets} -f GQ -Ob -o ${outdir}/snps.{}.bcf -" \
    ::: chr{1..14}
elif [[ ${mode} == 'mig' ]]; then
  # joint genotype calling
  parallel -j ${njobs} \
    "bcftools mpileup -b ${outdir}/bamlist -D -f ${ref} -q 30 -Q 30 -r {} -T ${dir_sites}/snps.sites.{} --ns UNMAP,SECONDARY,QCFAIL,DUP --lu PROPER_PAIR -a AD -Ou \
    | bcftools call -mv -G ${sample_sets} -f GQ -Ob -o ${outdir}/snps.{}.bcf -" \
    ::: chr{1..14}
else
  echo "Unrecognized argument: ${mode} ... expected 'phylo' or 'mig'"
  exit 1
fi

# concatenate chromosome BCFs
bcftools concat -Ob -o ${outdir}/snps.bcf ${outdir}/snps.chr{1..14}.bcf

# clean directory
rm ${outdir}/snps.chr{1..14}.bcf

# filter multisample BCF
bcftools filter -e "FMT/GQ<20" -S . -Ou ${outdir}/snps.bcf \
  | bcftools filter -e "DP<${min_dp} || DP>${max_dp} || MQ<30 || QUAL<30 || F_MISSING>0.1" -s 'LowQual' -Ou - \
  | bcftools view -m 2 -M 2 -v snps -c 1:minor -i 'FILTER="PASS"' -Ob -o ${outdir}/snps.filtered.bcf -

# index BCF
bcftools index ${outdir}/snps.filtered.bcf