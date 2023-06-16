#!/bin/bash
#
# Usage:
#   ld_decay_multi.sh <indir>
#
# Description:
#   Fit LD decay curves for multiple species with 'fit_LDdecay.R'.
#
# Requirements:
#   ngsld (incl. fit_LDdecay.R)
#   r

# path to 'fit_LDdecay.R'
fit_ld_decay=${HOME}/software/ngsLD/scripts/fit_LDdecay.R

# multispecies plot of LD decay
ls $1/*.ld.sampled \
  | awk 'BEGIN { print "File\tTaxon" } { sp=$1; sub(".*/", "", sp); sub("[.].*", "", sp); sub(/\w/, substr(toupper(sp), 1, 1), sp); print $1"\t"sp }' \
  | Rscript --vanilla --slave \
    ${fit_ld_decay} \
      --header \
      --ld r2 \
      --max_kb_dist 500 \
      --fit_boot 100 \
      --fit_bin_size 250 \
      --fit_level 100 \
      --plot_group 'Taxon' \
      --plot_scale 3 \
      --out $1/multispecies.ld_decay.pdf \
      &> $1/multispecies.fit_LDdecay.log