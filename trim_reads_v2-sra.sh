#!/bin/bash
#
# Usage:
#   trim_reads_v2-sra.sh <indir> <outdir_fastqs> <outdir_report> <njobs>
#
# Description:
#   Generates quality reports and perform trimming of FASTQ files with fastp.
#
# Requirements:
#   fastp
#   parallel
#
# Important:
#   Each job uses 8 CPU threads.

# trim reads and generate quality reports of before and after trimming
parallel -j $4 --link --plus \
  fastp \
    --in1 {1} \
    --in2 {2} \
    --out1 $2/{1/} \
    --out2 $2/{2/} \
    --detect_adapter_for_pe \
    --cut_right \
    --cut_right_window_size 4 \
    --cut_right_mean_quality 15 \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --n_base_limit 5 \
    --length_required 36 \
    --low_complexity_filter \
    --correction \
    --overrepresentation_analysis \
    --json $3/{1/..}-2.fastp.json \
    --html $3/{1/..}-2.fastp.html \
    --report_title {1/..}-2 \
    --thread 8 \
  ::: $1/{MA1,WOAK}*_1.fq.gz \
  ::: $1/{MA1,WOAK}*_2.fq.gz 