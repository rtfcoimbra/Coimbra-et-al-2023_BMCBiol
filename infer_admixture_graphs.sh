#!/bin/bash
#
# Usage:
#   infer_admixture_graphs.sh <treemix|orientagraph> <vcf> <plink.clst> <snp_block_size> <m_max> <n_rep> <outdir>
#
# Description:
#   Creates the input file for TreeMix from a VCF file and run either TreeMix or
#   OrientAGraph allowing for a range of migration edges with multiple replicates.
#
# Requirements:
#   orientagraph
#   plink
#   plink2treemix.py (https://bitbucket.org/nygcresearch/treemix/downloads/plink2treemix.py)
#   python2

mode=$1
vcf=$2
clst=$3
k=$4
m_max=$5
n_rep=$6
outdir=$7
out=${outdir}/$(basename ${vcf} .vcf)

# check if treemix input file exists
if [[ -f ${outdir}/treemix.frq.strat.gz ]]; then
  echo "INFO: TreeMix input file already exists!"
else
  # write a cluster-stratified frequency file
  plink \
    --vcf ${vcf} \
    --keep-allele-order \
    --chr-set 14 \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --freq gz \
    --missing \
    --within ${clst} \
    --out ${out}
  # convert stratified allele frequencies output from plink into treemix format
  python2 plink2treemix.py ${out}.frq.strat.gz ${outdir}/treemix.frq.strat.gz
fi

# check mode option
if [[ ${mode} != 'treemix' ]] && [[ ${mode} != 'orientagraph' ]]; then
  echo "Unrecognized argument: ${mode} ... expected 'treemix' or 'orientagraph'"
  exit 1
else
  mkdir -p ${outdir}/${mode}
  cd ${outdir}/${mode}
fi

# check running mode
if [[ ${mode} == 'treemix' ]]; then
  # run treemix allowing for a range of migration edges (m) with multiple replicates
  for m in $(seq ${m_max}); do
    for n in $(seq 1 ${n_rep}); do
      treemix \
        -i ${outdir}/treemix.frq.strat.gz \
        -o ${mode}.m${m}.r${n} \
        -root Outgroup \
        -k ${k} \
        -m ${m} \
        -bootstrap \
        -global \
        -se \
        -seed $RANDOM \
        > ${mode}.m${m}.r${n}.log
    done
  done
elif [[ ${mode} == 'orientagraph' ]]; then
  # run orientagraph allowing for a range of migration edges (m) with multiple replicates
  for m in $(seq ${m_max}); do
    for n in $(seq 1 ${n_rep}); do
      orientagraph \
        -i ${outdir}/treemix.frq.strat.gz \
        -o ${mode}.m${m}.r${n} \
        -root Outgroup \
        -k ${k} \
        -m ${m} \
        -bootstrap \
        -global \
        -se \
        -seed $RANDOM \
        -mlno \
        -allmigs \
        > ${mode}.m${m}.r${n}.log
    done
  done
fi

# find the highest likelihood runs
for m in $(seq ${m_max}); do
  logs=$(ls -v ${mode}.m${m}.r*.llik)
  best_run=$(( for log in ${logs}; do grep -Po 'Exiting.* \K(\-*\d+\.\d+) ' ${log}; done ) | sort -gr | head -1 | grep -f - ${logs} | grep -Po '\.r\K[0-9]+')
  echo ${mode}.m${m}.r${best_run} >> best_runs.list
done