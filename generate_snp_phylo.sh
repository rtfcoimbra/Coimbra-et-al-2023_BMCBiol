#!/bin/bash
#
# Usage:
#   generate_snp_phylo.sh <bcf> <outgroup_id> <outdir> <ncpus>
#
# Description:
#   Randomly subsample sites from a BCF file, convert VCF to PHYLIP, and
#   generate a SNP-based phylogeny.
#
# Requirements:
#   bcftools
#   iqtree
#   python3
#   vcflib
#   vcf2phylip.py (in the same directory as 'generate_snp_phylo.sh')

bcf=$1
prefix=$(basename ${bcf//filtered/sampled} .bcf)
outgroup=$2
outdir=$3
ncpus=$4

# randomly sample 1% of the SNP sites
bcftools view -Ov ${bcf} | vcfrandomsample -r 0.01 -p $RANDOM > ${outdir}/${prefix}.vcf

# convert VCF to PHYLIP
python3 vcf2phylip.py -i ${outdir}/${prefix}.vcf -o ${outgroup}

# generate PHYLIP without constant, partially constant, and ambiguously constant sites
iqtree -s ${outdir}/${prefix}.min4.phy -m MFP+ASC -T ${ncpus}

# perform model selection for models with ascertainment bias correction followed
# by tree inference with 1,000 replicates of the UFboot2 (with nearest neighbor
# interchange optimization) and the SH-aLRT
iqtree -s ${outdir}/${prefix}.min4.phy.varsites.phy -m MFP+ASC -B 1000 -alrt 1000 -bnni -T ${ncpus}