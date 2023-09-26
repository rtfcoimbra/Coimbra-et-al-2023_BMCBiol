# Code from: Genomic analysis reveals limited hybridization among three giraffe species in Kenya

Code used to analyze whole-genome sequencing data of giraffe in Coimbra *et al.* (2023):

- Coimbra RTF, Winter S, Muneza A, Fennessy S, Otiende M, Mijele D, Masiaine S, Stacy-Dawes J, Fennessy J, Janke A (2023) Genomic analysis reveals limited hybridization among three giraffe species in Kenya. *BMC Biology*, under review.

**Note:** The plotting scripts do not necessarily reproduce the figures exactly as shown in the paper. In some cases, I used a free image editing software, namely [Krita](https://krita.org/en/), to assemble independent plots and add or correct plot annotations.

## Workflow

- `workflow.txt`: describes the steps and the context in which the scripts described below were used for processing and analysing the whole-genome sequencing data of giraffe.

### Read quality control

- `trim_reads_v2.sh` and `trim_reads_v2-sra.sh`: trim paired-end reads with [fastp](https://github.com/OpenGene/fastp).

### Read mapping

- `map_reads_v2.sh`: map reads against the Masai giraffe assembly with [BWA](https://github.com/lh3/bwa) and sort the output BAMs with [samtools](https://github.com/samtools/samtools).
- `mark_duplicates.sh`: mark PCR/optical duplicate reads with [Picard MarkDuplicates](https://broadinstitute.github.io/picard/).
- `mapping_flagstats.sh`: caculate mapping statistics with [samtools](https://github.com/samtools/samtools).
- `realigner_target_creator.sh`: create list of target intervals for indel realignment with [GATK](https://software.broadinstitute.org/gatk/).
- `indel_realigner.sh`: perform local realignment around indels with [GATK](https://software.broadinstitute.org/gatk/).
- `clean_bams_v3.sh`: remove bad reads (flags 4, 256, 512, or 1024) from BAM files and keep only properly paired reads (flag 2) mapped to non-repetitive regions in autosomes with [samtools](https://github.com/samtools/samtools).

### Site depth statistics

- `site_depth_sample.sh`: calculate the mean site depth per sample with [mosdepth](https://github.com/brentp/mosdepth).
- `site_depth_global.sh`: calculate the global site depth for multiple giraffe BAMs with [sambamba](https://github.com/biod/sambamba).
- `site_depth_stats_v2.py`: calculate summary statistics from the global site depth distribution (i.e. 5th and 95th percentiles, median, and median absolute deviation) with Python (NumPy and SciPy).

### SNP and genotype calling

- `snp_calling_v2.sh`: estimate genotype likelihoods with [ANGSD](https://github.com/ANGSD/angsd).
- `genotype_calling.sh`: call genotypes with [bcftools](https://github.com/samtools/bcftools).

### Linkage pruning

- `calc_pairwise_ld.sh`: calculate pairwise linkage disequilibrium with [ngsLD](https://github.com/fgvieira/ngsLD).
- `ld_decay_multi.sh`: fit LD decay curves for multiple species with [fit_LDdecay.R](https://github.com/fgvieira/ngsLD/blob/master/scripts/fit_LDdecay.R).
- `ld_pruning.sh`: prune linked sites with [prune_graph.pl](https://github.com/fgvieira/ngsLD/blob/master/scripts/prune_graph.pl) and extract unlinked sites from ANGSD's genotype likelihoods file.
- `process_pruned_snps.sh`: combine pruned SNPs of all species, split them into separate files by autosome, and index resulting files for [ANGSD](https://github.com/ANGSD/angsd).

### Relatedness

- `infer_relatedness.sh`: infer pairwise relatedness among individuals based with [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) and [NGSremix](https://github.com/KHanghoj/NGSremix).
- `plot_figureS1.R`: plot relatedness coefficient as a heatmap and find the number of individuals to exclude (and identify them) for increasing thresholds of relatedness. ***Note:** this script was co-developed with Emma Vinson.*

### Population structure and admixture analyses

- `pcangsd_hwe.sh`: estimate covariance matrix and perform a Hardy-Weinberg equilibrium test with [PCAngsd](https://github.com/Rosemeis/pcangsd).
- `ngsadmix.sh`: estimate admixture proportions with [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix).
- `evaladmix.sh`: calculate a matrix of pairwise correlation of residuals between individuals using [evalAdmix](https://github.com/GenisGE/evalAdmix) to test the goodness of fit of the data to an admixture model.

### SNP-based phylogenomic inference

- `generate_snp_phylo.sh`: infer a phylogeny based on a random subset of SNPs from a BCF file using [bcftools](https://github.com/samtools/bcftools), [vcflib](https://github.com/vcflib/vcflib), [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip), and [IQ-TREE](http://www.iqtree.org/).

### Assembly and phylogeny of mitochondrial genomes

- See `workflow.txt`.

### Inference of migration events and test for introgression

- `infer_admixture_graphs.sh`: converts a VCF file into TreeMix input with [PLINK](https://www.cog-genomics.org/plink2/) and [`plink2treemix`](https://bitbucket.org/nygcresearch/treemix/downloads/plink2treemix.py), estimates admixture graphs with either [TreeMix](https://bitbucket.org/nygcresearch/treemix) or [OrientAGraph](https://github.com/sriramlab/OrientAGraph) allowing for a range of migration edges (*m*) with multiple bootstrap replicates, and finds the replicate with the highest likelihood per *m* value.
- `estimate_fbranch.sh`: calculate Patterson’s D, f4-ratio, and f-branch statistics with [Dsuite](https://github.com/millanek/Dsuite).

### Contemporary migration rates

- `run_ba3-snps.sh`: estimate comtemporary migration rates based on a random subset of SNPs from a BCF file using [bcftools](https://github.com/samtools/bcftools), [vcflib](https://github.com/vcflib/vcflib), [stacks](http://catchenlab.life.illinois.edu/stacks/), [stacksStr2immanc.pl](https://github.com/stevemussmann/file_converters/blob/master/stacksStr2immanc.pl), [BA3-SNPS-autotune](https://github.com/stevemussmann/BA3-SNPS-autotune) and [BA3-SNPS](https://github.com/stevemussmann/BayesAss3-SNPs).

### Demographic reconstruction

- `create_fasta.sh`: create genome consensus sequence in FASTA format with [sambamba](https://github.com/biod/sambamba) and [ANGSD](https://github.com/ANGSD/angsd).
- `estimate_sfs.sh`: estimate the unfolded site frequency spectrum (SFS) with [ANGSD](https://github.com/ANGSD/angsd) and realSFS.

## Figures

- **Figure 1:** [QGIS](https://www.qgis.org/en/site/) + `plot_figure1b-e.R` + `visFuns_modified.R` (script slightly modified from [`visFuns.R`](https://github.com/GenisGE/evalAdmix/blob/master/visFuns.R) from [evalAdmix](https://github.com/GenisGE/evalAdmix))
- **Figure 2:** `plot_figure2.R`
- **Figure 3:** `plot_figure3a.R` + `estimate_fbranch.sh` + `plot_figure3c.R`
- **Figure 4:** `plot_figure4.R`
- **Figure S1:** `plot_figureS1.R`
- **Figure S2:** `plot_figureS2.R`
- **Figure S3:** `plot_figureS3.R`
- **Figure S4:** `plot_figureS4.R` + `visFuns_modified.R` (script slightly modified from [`visFuns.R`](https://github.com/GenisGE/evalAdmix/blob/master/visFuns.R) from [evalAdmix](https://github.com/GenisGE/evalAdmix))
- **Figure S5:** `plot_figureS5.R`
- **Figure S6:** `plot_figureS6.R`
- **Figure S7:** `ld_decay_multi.sh`
- **Figure S8:** `plot_figureS8.R`