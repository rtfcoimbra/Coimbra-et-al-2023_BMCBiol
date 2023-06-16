################################################################################
#                                   SOFTWARE                                   #
################################################################################
#
# fastp             v0.20.0
# multiqc           v1.7
# seqkit            v0.11.0
# bwa               v0.7.17 (r1188)
# htslib            v1.10.2, v1.17
# samtools          v1.10
# picard            v2.21.7
# gatk              v3.8-1-0-gf15c1c3ef
# mosdepth          v0.3.3
# repeatmasker      v4.0.7-open
# bedtools          v2.27.1
# sambamba          v0.8.2
# angsd             v0.933, v0.940
# ngsld             v1.1.1
# ngsremix          v1.0.0
# pcangsd           v1.03
# ngsadmix          v33
# evaladmix         v0.962
# bcftools          v1.17
# vcflib            v1.0.3
# vcf2phylip        v2.8
# iqtree            v2.2.2.3
# getorganelle      v1.7.4
# mafft             v7.475
# plink             v1.9
# treemix           v1.13-r231
# optm              v0.1.6
# orientagraph      v1.0
# dsuite            v0.5-r52
# stacks            v2.41
# ba3-snps-autotune v2.1.2
# ba3-snps          v1.1
# stairway_plot     v2.1.1
# r                 v4.2.2
#
################################################################################
#                             DATA INTEGRITY CHECK                             #
################################################################################

# data storage directories
storage=/raw_data/mammals/Giraffa_spp/whole_genome_sequencing
sra_reads=/gendata_aj/rcoimbra/sra_reads

# check file content integrity
for dir in ${storage}/*; do
  cd ${dir}; md5sum --check MD5.txt >> ${storage}/md5sum_check.log
done

# create basic directory structure
mkdir \
  ~/project_kenya \
  ~/project_kenya/{data,results,scripts} \
  ~/project_kenya/data/{fasta,fastq,bam} \
  ~/project_kenya/data/fastq/{raw,trimmed}

################################################################################
#                      IDENTIFY REPEATS & PROCESS ASSEMBLY                     #
################################################################################

cd ~/project_kenya/data/fasta

# download chromosome-level Masai giraffe genome (Farré et al., 2019)
wget ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100590/giraffeChromosomes.fa.gz

# sort chromosomes by ID in natural order and add 'chr' to FASTA header
seqkit sort -N -n -i -2 giraffeChromosomes.fa.gz \
  | sed 's/>x/>X/g; s/>/>chr/g' > masai_giraffe_chromosomes.fa \
  && rm giraffeChromosomes.fa.gz

mkdir ~/project_kenya/results/repeatmasker

# identify repetitive regions in the reference genome assembly
RepeatMasker \
  -s \
  -engine 'ncbi' \
  -pa 32 \
  -species 'cetartiodactyla' \
  -dir ~/project_kenya/results/repeatmasker/ \
  -gff \
  ~/project_kenya/data/fasta/masai_giraffe_chromosomes.fa \
  > ~/project_kenya/results/repeatmasker/repeatmasker.log

# create a BED file of identified repeat regions while merging
# overlapping and adjacent repeats
tail -n +4 ~/project_kenya/results/repeatmasker/masai_giraffe_chromosomes.fa.out \
  | sed 's/^\s*//' \
  | sed -E 's/\s+/\t/g' \
  | cut -f 5-7 \
  | awk 'OFS="\t" { print $1, $2-1, $3 }' \
  | sort -V \
  | bedtools merge \
  > ~/project_kenya/results/repeatmasker/repeats.bed

# create a BED file of non-repetitive regions in autosomes
grep -v 'chrX' ~/project_kenya/data/fasta/masai_giraffe_chromosomes.fa.fai \
  | awk 'OFS="\t" { print $1,$2 }' \
  > ~/project_kenya/results/repeatmasker/autosomes.length
grep -v 'chrX' ~/project_kenya/results/repeatmasker/repeats.bed \
  | bedtools complement \
    -i stdin \
    -g ~/project_kenya/results/repeatmasker/autosomes.length \
  > ~/project_kenya/results/repeatmasker/autosomes.no_repeats.bed

# generate a new masked FASTA containing only autosomes
cut \
  -f 1 \
  ~/project_kenya/results/repeatmasker/autosomes.length \
  > ~/project_kenya/results/repeatmasker/autosomes.list
seqtk subseq \
  -l 80 \
  ~/project_kenya/results/repeatmasker/masai_giraffe_chromosomes.fa.masked \
  ~/project_kenya/results/repeatmasker/autosomes.list \
  > ~/project_kenya/results/repeatmasker/masai_giraffe_autosomes.masked.fa

# index new FASTA
samtools faidx ~/project_kenya/results/repeatmasker/masai_giraffe_autosomes.masked.fa

################################################################################
#                             READ QUALITY CONTROL                             #
################################################################################

# proceeed with data processing separately according to data source for
# easier handling of files

########################################################################
# process data from NovaSeq 6000 and HiSeq 4000

# set shell extglob
shopt -s extglob

# create symbolic links to raw FASTQs in local directory
ln -s ${storage}/!(ZNP01)/*.fq.gz ~/project_kenya/data/fastq/raw

# unset shell extglob
shopt -u extglob

mkdir ~/project_kenya/results/fastq_qc

# read quality control
~/project_kenya/scripts/trim_reads_v2.sh \
  ~/project_kenya/data/fastq/raw \
  ~/project_kenya/data/fastq/trimmed \
  ~/project_kenya/results/fastq_qc \
  7

# remove symlinks
for fq in ~/project_kenya/data/fastq/raw/*.fq.gz; do unlink ${fq}; done

########################################################################
# process data from HiSeq 2000

# create symbolic links to raw FASTQs in local directory
ln -s ${storage}/ZNP01/*.fq.gz ~/project_kenya/data/fastq/raw

# read quality control (add '--phred64')
r1=ZNP01_wHAIPI015412-108_C6L8JANXX_L5_1.fq.gz
r2=ZNP01_wHAIPI015412-108_C6L8JANXX_L5_2.fq.gz
fastp \
  --in1 ~/project_kenya/data/fastq/raw/${r1} \
  --in2 ~/project_kenya/data/fastq/raw/${r2} \
  --out1 ~/project_kenya/data/fastq/trimmed/${r1} \
  --out2 ~/project_kenya/data/fastq/trimmed/${r2} \
  --phred64 \
  --detect_adapter_for_pe \
  --cut_tail \
  --cut_tail_window_size 4 \
  --cut_tail_mean_quality 15 \
  --qualified_quality_phred 15 \
  --unqualified_percent_limit 40 \
  --n_base_limit 5 \
  --length_required 36 \
  --low_complexity_filter \
  --correction \
  --overrepresentation_analysis \
  --json ~/project_kenya/results/fastq_qc/${r1%.fq.gz}-2.fastp.json \
  --html ~/project_kenya/results/fastq_qc/${r1%.fq.gz}-2.fastp.html \
  --report_title ${r1%.fq.gz}-2 \
  --thread 8

# remove symlinks
for fq in ~/project_kenya/data/fastq/raw/*.fq.gz; do unlink ${fq}; done

########################################################################
# process data from SRA

# create symbolic links to raw FASTQs in local drectory
ln -s ${sra_reads}/*/*.fq.gz ~/project_kenya/data/fastq/raw

# read quality control (changed '--cut_tail' to '--cut_right')
~/project_kenya/scripts/trim_reads_v2-sra.sh \
  ~/project_kenya/data/fastq/raw \
  ~/project_kenya/data/fastq/trimmed \
  ~/project_kenya/results/fastq_qc \
  3

# remove symlinks
for fq in ~/project_kenya/data/fastq/raw/*.fq.gz; do unlink ${fq}; done

########################################################################

# generate read QC summary report
multiqc \
  --interactive \
  --filename 'multiqc_fastp.html' \
  --outdir ~/project_kenya/results/fastq_qc \
  ~/project_kenya/results/fastq_qc

################################################################################
#                                 READ MAPPING                                 #
################################################################################

cd ~/project_kenya/data/bam

mkdir ~/project_kenya/results/mapping

# read mapping
~/project_kenya/scripts/map_reads_v2.sh \
  ~/project_kenya/data/fasta/masai_giraffe_chromosomes.fa \
  ~/project_kenya/data/fastq/trimmed \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/mapping \
  40

# get names of samples sequenced in multiple lanes
ls *.sorted.bam | cut -d '_' -f 1 | sort | uniq -d > sample_ids.txt

# merge lane-level BAMs into sample-level BAMs
while read sample_id; do
  sample_bam=${sample_id}.sorted.bam
  echo "samtools merge ${sample_bam} ${sample_id}_*.sorted.bam" >> samtools-merge.jobs
done < sample_ids.txt
cat samtools-merge.jobs | parallel -j 8

# remove lane-level BAMs
while read sample_id; do
  rm ${sample_id}_*.sorted.bam
done < sample_ids.txt

# rename remaining BAMs
for bam in *_{L,SRR}*.sorted.bam; do
  mv ${bam} ${bam/_*.sorted.bam/.sorted.bam}
done

################################################################################
#                                MARK DUPLICATES                               #
################################################################################

mkdir ~/project_kenya/results/deduplication

# mark PCR/optical duplicate reads for unpatterned flowcells
java -XX:ParallelGCThreads=4 -Xmx20G -jar ~/software/picard.jar \
  MarkDuplicates \
    I=~/project_kenya/data/bam/ZNP01.sorted.bam \
    O=~/project_kenya/data/bam/ZNP01.dedup.bam \
    M=~/project_kenya/results/deduplication/ZNP01.dedup.metrics.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
    &> ~/project_kenya/results/deduplication/ZNP01.markduplicates.log

# remove sorted ZNP01 BAM
rm ~/project_kenya/data/bam/ZNP01.sorted.bam

# mark PCR/optical duplicate reads for patterned flowcells
~/project_kenya/scripts/mark_duplicates.sh \
  ~/project_kenya/data/bam \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/deduplication \
  5

# remove sorted BAMs
rm ~/project_kenya/data/bam/*.sorted.bam

################################################################################
#                             MAPPING QUALITY CHECK                            #
################################################################################

mkdir ~/project_kenya/results/mapping_flagstats

# caculate mapping statistics
~/project_kenya/scripts/mapping_flagstats.sh \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/mapping_flagstats \
  32

################################################################################
#                               INDEL REALIGNMENT                              #
################################################################################

mkdir ~/project_kenya/data/bam/{nubian,reticulated,masai,okapi}

# TODO CREATE LINKS FOR EACH SPECIES BAMS

mkdir -p ~/project_kenya/results/indel_realignment/{nubian,reticulated,masai,okapi}

# create separate lists of target intervals for each giraffe spp. and okapi
for spp in 'nubian' 'reticulated' 'masai' 'okapi'; do
  if [[ ${spp} == 'okapi' ]]; then
    mode='okapi'
  else
    mode='giraffe'
  fi
  ~/project_kenya/scripts/realigner_target_creator.sh \
    ${mode} \
    ~/project_kenya/data/fasta/masai_giraffe_chromosomes.fa \
    ~/project_kenya/data/bam/${spp} \
    ~/project_kenya/results/indel_realignment/${spp} \
    10 \
    4
done

# perform local realignment around indels
for spp in 'nubian' 'reticulated' 'masai' 'okapi'; do
  ~/project_kenya/scripts/indel_realigner.sh \
    ~/project_kenya/data/fasta/masai_giraffe_chromosomes.fa \
    ~/project_kenya/data/bam/${spp} \
    ~/project_kenya/results/indel_realignment/${spp} \
    ~/project_kenya/data/bam/${spp} \
    ~/project_kenya/results/indel_realignment/${spp} \
    5
done

# remove deduplicated BAMs
rm ~/project_kenya/data/bam/*.dedup.bam

# remove links to deduplicated BAMs
for spp in 'nubian' 'reticulated' 'masai' 'okapi'; do
  rm ~/project_kenya/data/bam/${spp}/*.dedup.bam
done

################################################################################
#                                  CLEAN BAMS                                  #
################################################################################

# remove reads mapped to repetitive regions and to sex chromosomes, and
# clean realigned BAMs
~/project_kenya/scripts/clean_bams_v3.sh \
  ~/project_kenya/results/repeatmasker/autosomes.no_repeats.bed \
  ~/project_kenya/data/bam \
  ~/project_kenya/data/bam \
  5

################################################################################
#                               DEPTH STATISTICS                               #
################################################################################

mkdir -p ~/project_kenya/results/depth_stats/{sample,global,nubian,reticulated,masai}

# calculate site depth per sample
~/project_kenya/scripts/site_depth_sample.sh \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/depth_stats/sample \
  20

# individuals excluded from GL- and GC-based analyses due to extremely
# low mean depth (< 1x): GF088 and GF184

# individuals excluded from GC-based analyses due to very low mean depth
# (< 6x): GF085, GF097, GF109, GF114, GF115, GF117, GF126, GF130, GF132,
# GF137, GF156, GF159, GF164, GF165, GF168, GF182, GF189, GF193, GF194,
# GF233, and ISC01

# iterate over species
for spp in 'masai' 'nubian' 'reticulated'; do
  # calculate global site depth
  ~/project_kenya/scripts/site_depth_global.sh \
    'gl' \
    ~/project_kenya/data/bam/${spp} \
    ~/project_kenya/results/depth_stats/${spp} \
    4 \
    5
  # calculate summary statistics for global site depth
  python3 ~/project_kenya/scripts/site_depth_stats_v2.py \
    ~/project_kenya/results/depth_stats/${spp}/site_depth.global.sampled > \
    ~/project_kenya/results/depth_stats/${spp}/site_depth.global.sampled.stats
done

# calculate global site depth
~/project_kenya/scripts/site_depth_global.sh \
  'gl' \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/depth_stats/global \
  4 \
  5

# calculate summary statistics for global site depth
python3 ~/project_kenya/scripts/site_depth_stats_v2.py \
  ~/project_kenya/results/depth_stats/global/site_depth.global.sampled > \
  ~/project_kenya/results/depth_stats/global/site_depth.global.sampled.stats

################################################################################
#                             SNP CALLING PER TAXON                            #
################################################################################

mkdir -p ~/project_kenya/results/snp_calling/{masai,nubian,reticulated}

# iterate over species
for spp in 'masai' 'nubian' 'reticulated'; do
  # call and filter SNPs
  ~/project_kenya/scripts/snp_calling_v2.sh \
    ~/project_kenya/results/repeatmasker/masai_giraffe_autosomes.masked.fa \
    ~/project_kenya/data/bam/${spp} \
    ~/project_kenya/results/depth_stats/${spp}/site_depth.global.sampled.stats \
    ~/project_kenya/results/snp_calling/${spp} \
    5
done

################################################################################
#                            LINKAGE DISEQUILIBRIUM                            #
################################################################################

mkdir -p ~/project_kenya/results/ld_pruning/{masai,nubian,reticulated}

# calculate pairwise LD per species
for spp in 'masai' 'nubian' 'reticulated'; do
  ~/project_kenya/scripts/calc_pairwise_ld.sh \
    ~/project_kenya/results/snp_calling/${spp} \
    ~/project_kenya/results/ld_pruning/${spp} \
    32
done

# fit LD decay curves
for spp in 'masai' 'nubian' 'reticulated'; do
  mv \
    ~/project_kenya/results/ld_pruning/${spp}/snps.ld.sampled \
    ~/project_kenya/results/ld_pruning/${spp}.snps.ld.sampled
done
~/project_kenya/scripts/ld_decay_multi.sh ~/project_kenya/results/ld_pruning

# iterate over species
for spp in 'masai' 'nubian' 'reticulated'; do
  # set appropriate minimum weight for LD pruning
  if [[ ${spp} == 'nubian' ]]; then
    min_weight=0.15
  else
    min_weight=0.1
  fi
  # prune SNPs for linkage disequilibrium
  ~/project_kenya/scripts/ld_pruning.sh \
    ~/project_kenya/results/snp_calling/${spp} \
    ~/project_kenya/results/ld_pruning/${spp} \
    100 \
    ${min_weight} \
    7
done

# combine pruned SNPs of all species and prepare sites files per
# autosome for angsd
~/project_kenya/scripts/process_pruned_snps.sh \
  ~/project_kenya/results/ld_pruning \
  14

################################################################################
#                             SNP CALLING COMBINED                             #
################################################################################

# differences compared to SNP calling per taxon: site depth filters are 
# based on estimates of combined global depth; no filters for minMAF,
# sb_pval, hetbias_pval, hwe_pval, and snp_pval are used.

mkdir ~/project_kenya/results/snp_calling/all_pruned

~/project_kenya/scripts/snp_calling_v2.sh \
  'all' \
  ~/project_kenya/results/repeatmasker/masai_giraffe_autosomes.masked.fa \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/depth_stats/global/site_depth.global.sampled.stats \
  ~/project_kenya/results/ld_pruning \
  ~/project_kenya/results/snp_calling/all_pruned \
  5

################################################################################
#                             FILTER BY RELATEDNESS                            #
################################################################################

mkdir ~/project_kenya/results/relatedness

# estimate individual ancestries, allele frequencies, and relatedness assuming K=3
~/project_kenya/scripts/infer_relatedness.sh \
  ~/software/maf_beagle_ngsadmix/maf_beagle \
  ~/project_kenya/results/snp_calling/all_pruned \
  ~/project_kenya/results/relatedness \
  3 \
  0.001 \
  10

# run script plot_figureS1.R

# remove individuals: GF096, GF132, GF209, GF217, GF226, GF231, GF233, GF257,
# GF259, and GF262

################################################################################
#                        DEPTH STATISTICS FOR UNRELATED                        #
################################################################################

mkdir ~/project_kenya/results/depth_stats/{global_unrelated,nubian_unrelated,reticulated_unrelated,masai_unrelated}

# iterate over species
for spp in 'masai' 'nubian' 'reticulated'; do
  # calculate global site depth
  ~/project_kenya/scripts/site_depth_global.sh \
    'gl' \
    ~/project_kenya/data/bam/${spp} \
    ~/project_kenya/results/depth_stats/${spp}_unrelated \
    4 \
    5
  # calculate summary statistics for global site depth
  python3 ~/project_kenya/scripts/site_depth_stats_v2.py \
    ~/project_kenya/results/depth_stats/${spp}_unrelated/site_depth.global.sampled > \
    ~/project_kenya/results/depth_stats/${spp}_unrelated/site_depth.global.sampled.stats
done

# calculate global site depth for unrelated individuals
~/project_kenya/scripts/site_depth_global.sh \
  'gl' \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/depth_stats/global_unrelated \
  4 \
  5

# calculate summary statistics for global site depth of unrelated individuals
python3 ~/project_kenya/scripts/site_depth_stats_v2.py \
  ~/project_kenya/results/depth_stats/global_unrelated/site_depth.global.sampled > \
  ~/project_kenya/results/depth_stats/global_unrelated/site_depth.global.sampled.stats

################################################################################
#                      SNP CALLING COMBINED FOR UNRELATED                      #
################################################################################

mkdir ~/project_kenya/results/snp_calling/unrelated_pruned

# joint snp calling for unrelated individuals
~/project_kenya/scripts/snp_calling_v2.sh \
  'all' \
  ~/project_kenya/results/repeatmasker/masai_giraffe_autosomes.masked.fa \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/depth_stats/global_unrelated/site_depth.global.sampled.stats \
  ~/project_kenya/results/ld_pruning \
  ~/project_kenya/results/snp_calling/unrelated_pruned \
  5

################################################################################
#                             POPULATION STRUCTURE                             #
################################################################################

mkdir ~/project_kenya/results/{pca,admixture,evaladmix}

# calculate a covariance matrix for PCA
~/project_kenya/scripts/pcangsd_hwe.sh \
  ~/project_kenya/results/snp_calling/unrelated_pruned \
  ~/project_kenya/results/pca \
  0.001 \
  10

# estimate admixture proportions
~/project_kenya/scripts/ngsadmix.sh \
  ~/project_kenya/results/snp_calling/unrelated_pruned \
  ~/project_kenya/results/admixture \
  1 \
  14 \
  0.001 \
  20

# assess model fit of admixture models
~/project_kenya/scripts/evaladmix.sh \
  ~/project_kenya/results/snp_calling/unrelated_pruned \
  ~/project_kenya/results/admixture \
  ~/project_kenya/results/evaladmix \
  1 \
  14 \
  0.001 \
  10

################################################################################
#                                 SNP PHYLOGENY                                #
################################################################################

mkdir ~/project_kenya/results/depth_stats/{global_unrelated_snp-phylo}

# calculate global site depth for all giraffe and okapi
~/project_kenya/scripts/site_depth_global.sh \
  'phylo' \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/depth_stats/global_unrelated_snp-phylo \
  4 \
  5

# calculate summary statistics for global site depth of all giraffe and okapi
python3 ~/project_kenya/scripts/site_depth_stats_v2.py \
  ~/project_kenya/results/depth_stats/global_unrelated_snp-phylo/site_depth.global.sampled > \
  ~/project_kenya/results/depth_stats/global_unrelated_snp-phylo/site_depth.global.sampled.stats

mkdir ~/project_kenya/results/genotype_calling_snp-phylo

# call and filter genotypes
~/project_kenya/scripts/genotype_calling.sh \
  'phylo' \
  ~/project_kenya/results/repeatmasker/masai_giraffe_autosomes.masked.fa \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/depth_stats/global_unrelated_snp-phylo/site_depth.global.stats \
  ~/project_kenya/results/genotype_calling_snp-phylo \
  7

mkdir ~/project_kenya/results/snp_phylogeny

# randomly sample 1% of all SNPs and generate SNP-phylogeny
~/project_kenya/scripts/generate_snp_phylo.sh \
  ~/project_kenya/results/genotype_calling_snp-phylo/snps.filtered.bcf \
  WOAK \
  ~/project_kenya/results/snp_phylogeny \
  20

################################################################################
#                         INFERENCE OF MIGRATION EVENTS                        #
################################################################################

mkdir ~/project_kenya/results/admixture_graphs

# run treemix to obtain likelihoods for a broad range of allowed migration edges (m)
~/project_kenya/scripts/infer_admixture_graphs.sh \
  'treemix' \
  ~/project_kenya/results/snp_phylogeny/snps.sampled.vcf \
  ~/project_kenya/results/admixture_graphs/sets.clst \
  100 \
  5 \
  50 \
  ~/project_kenya/results/admixture_graphs

# run script plot_figureS5.R to narrow down m values for further investigation

# run orientagraph to obtain reliable graphs for the narrowed down range of m values
~/project_kenya/scripts/infer_admixture_graphs.sh \
  'orientagraph' \
  ~/project_kenya/results/snp_phylogeny/snps.sampled.vcf \
  ~/project_kenya/results/admixture_graphs/sets.clst \
  100 \
  2 \
  10 \
  ~/project_kenya/results/admixture_graphs

################################################################################
#                         INFERENCE OF MIGRATION EVENTS                        #
################################################################################

# dsuite


################################################################################
#                            RECENT MIGRATION RATES                            #
################################################################################

mkdir ~/project_kenya/results/depth_stats/{global_unrelated_ba3-snps}

# calculate global site depth for unrelated individuals
~/project_kenya/scripts/site_depth_global.sh \
  'mig' \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/depth_stats/global_unrelated_ba3-snps \
  4 \
  5

# calculate summary statistics for global site depth of unrelated individuals
python3 ~/project_kenya/scripts/site_depth_stats_v2.py \
  ~/project_kenya/results/depth_stats/global_unrelated_ba3-snps/site_depth.global.sampled > \
  ~/project_kenya/results/depth_stats/global_unrelated_ba3-snps/site_depth.global.sampled.stats

mkdir ~/project_kenya/results/genotype_calling_ba3-snps

# call and filter genotypes
~/project_kenya/scripts/genotype_calling.sh \
  'mig' \
  ~/project_kenya/results/repeatmasker/masai_giraffe_autosomes.masked.fa \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/depth_stats/global_unrelated_ba3-snps/site_depth.global.stats \
  ~/project_kenya/results/genotype_calling_ba3-snps \
  7

mkdir ~/project_kenya/results/recent_gene_flow

# infer recent migration rates
~/project_kenya/scripts/run_ba3-snps.sh \
  ~/project_kenya/results/genotype_calling_ba3-snps/snps.filtered.bcf \
  ~/project_kenya/results/genotype_calling_ba3-snps/sample_sets.txt \
  0.02 \
  3 \
  22000000 \
  2000000 \
  2000 \
  ~/project_kenya/results/recent_gene_flow \
  3

################################################################################
#                              RECENT DEMOGRAPHY                               #
################################################################################

mkdir -p ~/project_kenya/results/demography/{fasta,sfs,stairway_plot}

# generate okapi FASTA for polarizing the SFS
~/project_kenya/scripts/create_fasta.sh \
  ~/project_kenya/results/repeatmasker/masai_giraffe_autosomes.masked.fa \
  ~/project_kenya/data/bam \
  ~/project_kenya/results/demography/fasta \
  8

# estimate SFS per species
for spp in 'masai' 'nubian' 'reticulated'; do
  ~/project_kenya/scripts/estimate_sfs.sh \
    ~/project_kenya/results/repeatmasker/masai_giraffe_autosomes.masked.fa \
    ~/project_kenya/results/demography/fasta/WOAK.fa \
    ~/project_kenya/data/bam/${spp} \
    ~/project_kenya/results/depth_stats/${spp}_unrelated/site_depth.global.sampled.stats \
    ~/project_kenya/results/demography/sfs/${spp} \
    16
done

# run stairway plot without singletons