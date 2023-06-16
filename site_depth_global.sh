#!/bin/bash
#
# Usage:
#   site_depth_global.sh <phylo|mig|gl> <dir_bams> <outdir> <ncpus> <njobs>
#
# Description:
#   Calculate the global site depth from cleaned BAMs using sambamba.
#
# Requirements:
#   parallel
#   pigz
#   sambamba

mode=$1
dir_bams=$2
outdir=$3
ncpus=$4
njobs=$5

# set shell extglob
shopt -s extglob

# check running mode
if [[ ${mode} == 'phylo' ]]; then
  # create an array of all giraffe species and okapi BAMs for GC
  bams=(${dir_bams}/{GNP01,GNP04,GNP05,SNR2,ZNP01,ETH1,ETH2,ETH3,GF084,GF089,GF090,GF094,GF101,GF105,MF06,MF22,MF24,GF220,GF221,GF224,GF258,GF260,GF261,GF263,GF264,GF267,GF268,GF269,GF270,GF271,GF272,GF081,GF095,WA720,WA733,WA746,WA806,WA808,ENP11,ENP16,ENP19,ENP20,HNB102,HNB110,BNP02,KKR01,KKR08,MTNP09,V23,SUN3,GF289,GF211,GF213,GF214,GF290,GF227,GF228,GF229,GF230,GF232,GF234,GF235,GF236,GF237,ISC04,ISC08,GF292,GF295,GF276,GF277,GF278,GF281,GF283,GF286,GF288,GF187,GF201,GF203,GF206,GF218,GF223,LVNP8-04,LVNP8-08,LVNP8-09,LVNP8-10,LVNP8-12,LVNP8-14,GF155,GF080,GF048,GF054,GF070,GF074,GF076,GF250,GF253,GF068,GF146,GF242,GF245,GF248,GF249,GF049,GF077,GF246,MA1,GF005,GF007,GF011,GF016,GF023,GF034,GF035,GF036,GF238,GF239,GF044,GF065,GF144,SGR01,SGR05,SGR07,SGR13,SGR14,GF059,WOAK}.clean.bam)
elif [[ ${mode} == 'mig' ]]; then
  # create an array of nubian, reticulated, and masai giraffe BAMs for GC
  bams=(${dir_bams}/{nubian,reticulated,masai}/!(GF085|GF097|GF109|GF114|GF115|GF117|GF126|GF130|GF132|GF137|GF156|GF159|GF164|GF165|GF168|GF182|GF189|GF193|GF194|GF233|ISC01).clean.bam)
elif [[ ${mode} == 'gl' ]]; then
  # create an array of nubian, reticulated, and masai giraffe BAMs for GL
  bams=(${dir_bams}/!(WOAK).clean.bam)
else
  echo "Unrecognized argument: ${mode} ... expected 'phylo', 'mig' or 'gl'"
  exit 1
fi

# unset shell extglob
shopt -u extglob

# calculate global site depth per autosome in parallel
parallel -j ${njobs} \
  "sambamba depth base -L {} -t ${ncpus} --combined --fix-mate-overlaps ${bams[@]} | awk 'NR>1 { print \$3 }' > ${outdir}/{}.site_depth.global" \
  ::: chr{1..14}

# concatenate autosomes site depth estimates into a single file
cat ${outdir}/chr{1..14}.site_depth.global > ${outdir}/site_depth.global

# randomly sample 0.01% of all sites for plotting the global site depth distribution
awk 'rand()<0.0001' ${outdir}/site_depth.global > ${outdir}/site_depth.global.sampled

# compress files
pigz ${outdir}/*.global