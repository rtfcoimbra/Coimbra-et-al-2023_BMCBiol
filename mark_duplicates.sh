#!/bin/bash
#
# Usage:
#   mark_duplicates.sh <indir> <outdir_bam> <outdir_log> <njobs>
#
# Description:
#   Mark PCR/optical duplicate reads with Picard's MarkDuplicates.
#
# Requirements:
#   parallel
#   picard
#   samtools
#
# Important:
#   Each Picard job uses 4 CPU threads.

# path to Picard
picard=~/software/picard.jar

# mark PCR/optical duplicate reads for patterned flowcells
parallel -j $4 --plus \
  "java -XX:ParallelGCThreads=4 -Xmx20G -jar ${picard} \
    MarkDuplicates \
      I={} \
      O=$2/{/..}.dedup.bam \
      M=$3/{/..}.dedup.metrics.txt \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      &> $3/{/..}.markduplicates.log" \
  ::: $1/*.sorted.bam

# index new BAMs
parallel -j $4 samtools index -b {} ::: $2/*.dedup.bam
