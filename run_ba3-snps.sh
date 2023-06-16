#!/bin/bash
#
# Usage:
#   run_ba3-snps.sh <bcf> <popmap> <snp_sampling_prob> <nruns> <niter> <burnin> <sampling> <outdir> <njobs>
#
# Description:
#   Randomly subsample sites from a BCF file, convert VCF to IMMANC, autotune
#   mixing parameters for BA3-SNP, and estimate comtemporary migration rates
#   with multiple parallel runs of BA3-SNPs.
#
# Requirements:
#   BA3-SNPS
#   BA3-SNPS-autotune
#   bcftools
#   parallel
#   perl
#   stacks
#   stacksStr2immanc.pl
#   vcflib

bcf=$1
popmap=$2
prob=$3
nruns=$4
niter=$5
burnin=$6
sampling=$7
outdir=$8
njobs=$9

prefix=${outdir}/$(basename ${bcf//filtered/sampled} .bcf)
vcf=${prefix}.vcf
immanc=${prefix}.immanc

# randomly sample a fraction of all SNP sites
bcftools view -Ov ${bcf} | vcfrandomsample -r ${prob} > ${vcf}

# convert VCF to IMMANC format
populations -V ${vcf} -O ${outdir} -M ${popmap} --structure
perl stacksStr2immanc.pl -s ${prefix}.p.structure -o ${immanc}

# check missing data frequency per individual
start='BEGIN variant_sites_per_sample'; end='END variant_sites_per_sample'
sed -n "/${start}/, /${end}/{ /${start}/! { /${end}/! p } }" ${prefix}.p.log.distribs \
  | tail -n +3 \
  | cut -f 1,5 \
  | sort -grk 2 > ${prefix}.missfreq

# autotune mixing parameters for BA3-SNPS (i.e. M, A, F)
nsites=$(bcftools view -H ${vcf} | wc -l)
BA3-SNPS-autotune.py -i ${immanc} -l ${nsites} -b 500000 -g 2500000 -o ${prefix}
read -r M A F <<< $(tail -n +3 ${immanc}.finalParams.txt)

# organize files
mkdir ${outdir}/{autotune,stacks}
mv ${prefix}.p.* ${outdir}/stacks
mv ${prefix}{,*.stdout,*.txt} ${outdir}/autotune

# run BA3-SNPS
for i in $(seq ${nruns}); do
  immanc_ln=${outdir}/run${i}/$(basename ${immanc})
  out=${immanc_ln%.immanc}.run${i}
  echo "mkdir ${outdir}/run${i}; cd ${outdir}/run${i}; ln -s ${immanc} ${immanc_ln}; BA3-SNPS -F ${immanc_ln} -l ${nsites} -s $RANDOM -i ${niter} -b ${burnin} -n ${sampling} -m ${M} -a ${A} -f ${F} -v -u -t -o ${out} > ${out}.log" >> ${outdir}/ba3-snps.jobs
done
cat ${outdir}/ba3-snps.jobs | parallel -j ${njobs}