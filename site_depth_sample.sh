#!/bin/bash
#
# Usage:
#   site_depth_sample.sh <indir> <outdir> <njobs>
#
# Description:
#   Calculate per-base depth, and genome-wide and chromosome-wide cumulative
#   depth distributions from clean BAMs using mosdepth. Reads mapped in a
#   proper pair (2) are included in calculations, while unmapped (4),
#   secondary (256), QC failed (512), and duplicate (1024) reads are excluded.
#
# Requirements:
#   mosdepth
#   parallel
#
# Important:
#   Each job uses 4 CPU threads.

# calculate the cumulative depth distribution
parallel -j $3 --plus \
  mosdepth -t 4 -F 1796 -i 2 $2/{/..} {} \
  ::: $1/*.clean.bam
