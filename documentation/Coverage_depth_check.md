# Coverage depth check for somatic variant calling	

## Overview

To analyse the evenness of coverage depth of the samples, we used a file that contained the coordinates of the baits used for exome capture. Coordinates from hg19 were lifovered to GRCh38. Following this, coverage observed was performed using `samtools depth` command to calculate the coverage depth of each sample. The output was then processed to generate a summary of the coverage depth across all samples.

Samples with a coverage depth of **>20X across at least 80% of the baits** used were considered for somatic variant calling. 

## Dependencies

## Required software

The scripts used for this analysis are located in the `scripts/coverage_qc` directory. The analysis was performed using a combination of R, Perl, and shell scripts. The following software versions were used:

- Perl scripts were run using Perl version `v5.38.0`[**here**](https://www.perl.org/)
- `LifOver` [**here**](https://genome.ucsc.edu/cgi-bin/hgLiftOver)
- `samtools` version `1.14` [**here**](https://www.htslib.org/)
- `R` version `4.3.3` [**here**](https://www.r-project.org/)

### Required datasets

- The bait set used for the exome capture was the Agilent SureSelect Human All Exon V6+UTR, the file was obtained from Agilent (file S07604624_hs_hg19.zip). Coordinates for GRCh38 were obtained by performing a liftover from hg19 to GRCh38 coordinates using [`LifOver`](https://genome.ucsc.edu/cgi-bin/hgLiftOver). See file [`resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed`](resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed) 

## Coverage Depth analysis

The Coverage Depth analysis was run under internal cluster environment with **Ubuntu 22.04** and **LSF job execution**. It is required to replace the `PROJECTDIR` shell environment variable with the path of where the current repository was downloaded. 

- To make a copy Links to the bam files I ran the following 
```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
cp -L /nfs/cancer_ref01/nst_links/live/2744/*/*.sample.dupmarked.bam ${PROJECTDIR}/BAMS/
#Then to Copy the indexes I ran the following
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
cp -L /nfs/cancer_ref01/nst_links/live/2744/*/*.sample.dupmarked.bam.bai ${PROJECTDIR}/BAMS/
```
### To make a directory to store the results of the coverage depth analysis

- To make a job array to calculate the depth with Samtools 

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes

mkdir -p ${PROJECTDIR}/results/qc_plots/depth
mkdir -p ${PROJECTDIR}/logs/depth
cd ${PROJECTDIR}/results/qc_plots/depth
#To create the sample list file 
ls -1 ${PROJECTDIR}/BAMS/*.bam >sample_bams.fofn
ls -1 ${PROJECTDIR}/BAMS/*.bam | cut -f 9 -d "/" | sed 's/\.sample\.dupmarked\.bam//' >sample.list
#To create the jobs for submission of sample depth calculation
for f in `cat sample.list`; do bsub -e ${PROJECTDIR}/logs/depth/depth.$f.e -o ${PROJECTDIR}/logs/depth/depth.$f.o -n 2 -M2000 -R"select[mem>2000] rusage[mem=2000]" "module load samtools/1.14; samtools depth -a -b ${PROJECTDIR}/required_datasets/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed -o $f.depth.tsv -J -s -@ 2 ${PROJECTDIR}/BAMS/$f.sample.dupmarked.bam"; done

```
- Then to  summarise the coverage stats per sample

```bash
mkdir -p ${PROJECTDIR}/scripts
cp /lustre/scratch119/casm/team113da/projects/6500_DERMATLAS_SG_basal_cell_adenoma_and_adenocarcinoma_WES/SCRIPTS/DEPTH/* ${PROJECTDIR}/scripts
cd ${PROJECTDIR}/qc_plots/depth
for f in `cat sample.list`; do bsub -e ${PROJECTDIR}/logs/depth/count.$f.e -o ${PROJECTDIR}/logs/depth/count.$f.o -M2000 -R"select[mem>2000] rusage[mem=2000]" -q small "${PROJECTDIR}/scripts/count_region_coverage.pl $f.depth.tsv > $f.covstats.tsv"; done
```
-Then to summarise the coverage stats per cohort
```bash
#Togenerate the cov_stats summary for the cohort

#Get the header from the first file
for f in `ls -1 *.covstats.tsv | head -n1`; do head -n1 $f | sed 's/^11/Sample\t11/' >cov_stats_summary.tsv; done
grep -v Mean *stats.tsv | sed 's/.covstats.tsv:/\t/' >> cov_stats_summary.tsv

Rscript scripts/coverage_qc/02_plot_cov_stats_summary_sortByMeanCov_mod.R
```
