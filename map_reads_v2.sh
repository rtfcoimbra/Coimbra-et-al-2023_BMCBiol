#!/bin/bash
#
# Usage:
#   map_reads_v2.sh <ref.fa> <indir> <outdir_bam> <outdir_log> <ncpus>
#
# Description:
#   Map reads against a reference assembly using BWA-MEM, fill in mate
#   coordinates, insert size and mate related flags, and sort the output
#   BAMs with samtools.
#
# Requirements:
#   bwa
#   samtools

# index assembly
bwa index -a bwtsw $1

# iterate over read 1 files
for r1 in $2/*_1.fq.gz; do
  # read 2
  r2=${r1/_1./_2.}
  # basename
  base=$(basename ${r1})
  # output sorted BAM
  sorted_bam=$3/${base/_1.fq.gz/.sorted.bam}
  # logfile
  log=$4/${base/_1.fq.gz/.mapping.log}
  # generate read group fields
  if [[ ${base} == @(MA1|WOAK)* ]]; then #https://stackoverflow.com/questions/35812293/what-is-syntax-in-bash
    id=$(echo ${base} | awk -F '_' '{ print $1"."$2 }')
    sm=$(echo ${base} | awk -F '_' '{ print $1 }')
    lb=$(echo ${base} | awk -F '_' '{ print $1".SRR" }')
  else
    id=$(echo ${base} | sed -e 's/_/\t/g; s/L\([0-9]\)/\1/g' | awk '{ print $1"."$3"."$4 }')
    sm=$(echo ${base} | awk -F '_' '{ print $1 }')
    lb=$(echo ${base} | awk -F '_' '{ print $2 }')
  fi
  # map reads against reference; fill in mate coordinates, insert size and mate
  # related flags; sort BAM by coordinates
  { bwa mem -M -t $5 -R "@RG\\tID:${id}\\tSM:${sm}\\tPL:ILLUMINA\\tLB:${lb}" $1 ${r1} ${r2} \
    | samtools fixmate -@ 4 - - \
    | samtools sort -@ 4 -o ${sorted_bam} - ; } 2> ${log}
done
