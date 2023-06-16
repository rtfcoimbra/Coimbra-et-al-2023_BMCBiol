#!/bin/bash
#
# Usage:
#   create_fasta.sh <ref.fa> <bam> <outdir> <njobs>
#
# Description:
#   Calculate per-base depth from clean BAMs using sambamba and create a
#   consensus sequence in FASTA format with ANGSD.
#
# Requirements:
#   angsd
#   parallel
#   python3
#   sambamba
#   samtools
#   'site_depth_stats_v2.py' (must be in same directory as 'generate_fasta.sh')
#
# Important:
#   <ref.fa> must be the masked Masai giraffe assembly with sex chromosomes
#   removed, i.e. the 'masai_giraffe_autosomes.masked.fa'.

ref=$1
bam=$2
sample=$(basename ${bam%.clean.bam})
outdir=$3
njobs=$4

# calculate site depth per autosome in parallel
parallel -j ${njobs} \
  "sambamba depth base -L {} -t 4 --fix-mate-overlaps ${bam} | awk 'NR>1 { print \$3 }' > ${outdir}/{}.site_depth" \
  ::: chr{1..14}

# concatenate autosomes site depth estimates into a single file
cat ${outdir}/chr{1..14}.site_depth > ${outdir}/${sample}.site_depth
&& rm ${outdir}/chr{1..14}.site_depth

# calculate site depth statistics
python3 site_depth_stats_v2.py ${outdir}/${sample}.site_depth ${outdir}/${sample}.site_depth.stats

# get the 95th percentile of the site depth distribution
max_dp=$(grep -Po '95th percentile: \K\d+' ${outdir}/${sample}.site_depth.stats)
# set read filters
filter_reads="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -baq 2"
# set site filters
filter_sites="-minMapQ 30 -minQ 30 -setMinDepthInd 4 -setMaxDepthInd ${max_dp}"
# set angsd tasks
todo="-doCounts 1 -doFasta 1 -basesPerLine 80"

# generate FASTA file
angsd \
  -i ${bam} \
  -ref ${ref} \
  ${filter_reads} \
  ${filter_sites} \
  ${todo} \
  -P 4 \
  -out ${outdir}/${sample} \
  &> ${outdir}/${sample}.log

# decompress FASTA
gunzip ${outdir}/${sample}.fa.gz

# index FASTA
samtools faidx ${outdir}/${sample}.fa