#!/bin/bash
#
# Usage:
#   run_dsuite.sh <vcf> <tree.nwk> <sample_sets.txt> <outdir>
#
# Description:
#   
#
# Requirements:
#   dsuite
#   plot_d.rb
#   plot_f4ratio.rb
#   python3
#   ruby

Dsuite Dtrios -c -g -k 100 -o dsuite -n oag-m1 -t oag.m1.r2.rooted.nwk snps.sampled.vcf sample_sets.txt

#ruby plot_d.rb dsuite_oag-m1_BBAA.txt plot_order.txt 0.25 dsuite_oag-m1_BBAA_D.svg

#ruby plot_f4ratio.rb dsuite_oag-m1_BBAA.txt plot_order.txt 0.5 dsuite_oag-m1_BBAA_f4ratio.svg

Dsuite Fbranch oag.m1.r2.rooted.nwk dsuite_oag-m1_tree.txt > dsuite_oag-m1_Fbranch.txt

python3 ~/software/Dsuite/utils/dtools.py --dpi 300 dsuite_oag-m1_Fbranch.txt oag.m1.r2.rooted.nwk