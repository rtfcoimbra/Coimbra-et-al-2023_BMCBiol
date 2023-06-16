#!/bin/bash
#
# Usage:
#   mapping_flagstats.sh <indir> <outdir> <ncpus>
#
# Description:
#   Calculate mapping statistics with samtools.
#
# Requirements:
#   parallel
#   samtools

# calculate mapping statistics
parallel -j $3 "samtools flagstat {} > $2/{/.}.flagstats" ::: $1/*.dedup.bam