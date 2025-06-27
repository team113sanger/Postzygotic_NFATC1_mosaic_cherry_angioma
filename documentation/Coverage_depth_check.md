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

### Download BAM files 

- BAM files are available for Download from EGA repository [**here**](https://ega-archive.org/datasets/EGAD00001008664). 

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes

```
### Calculate depth of coverage

To make a job array to calculate the depth of coverage for each sample across the baitset coordinates we ran the following commands:

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
RESULTSDIR=${PROJECTDIR}/results/qc_plots/depth
BAITSET=${PROJECTDIR}/resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed

mkdir -p ${RESULTSDIR}
mkdir -p ${PROJECTDIR}/logs/depth

cd ${RESULTSDIR}
#To create the sample list file 
ls -1 ${PROJECTDIR}/BAMS/*.bam >sample_bams.fofn
ls -1 ${PROJECTDIR}/BAMS/*.bam | cut -f 9 -d "/" | sed 's/\.sample\.dupmarked\.bam//' >sample.list
#To create the jobs for submission of sample depth calculation
for f in `cat sample.list`; do bsub -e ${PROJECTDIR}/logs/depth/depth.${f}.e -o ${PROJECTDIR}/logs/depth/depth.${f}.o -n 2 -M2000 -R"select[mem>2000] rusage[mem=2000]" "samtools depth -a -b ${BAITSET} -o ${f}.depth.tsv -J -s -@ 2 ${PROJECTDIR}/BAMS/${f}.sample.dupmarked.bam"; done
```
OUTPUT : `${f}.depth.tsv` files for each sample containing 

Subsequently to generate the summaries of the proportion of bases covered with at least a given coverage threshold in steps of 10X, the following commands were run:

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
RESULTSDIR=${PROJECTDIR}/results/qc_plots/depth
DPSCRIPTDIR=${PROJECTDIR}/scripts/coverage_qc

cd ${RESULTSDIR}

for f in `cat sample.list`; do bsub -e ${PROJECTDIR}/logs/depth/count.${f}.e -o ${PROJECTDIR}/logs/depth/count.${f}.o -M2000 -R"select[mem>2000] rusage[mem=2000]" -q small "${DPSCRIPTDIR}/count_region_coverage.pl ${f}.depth.tsv > ${f}.covstats.tsv"; done
```
OUTPUT : `${f}.covstats.tsv` files for each sample containing the coverage depth statistics for each sample at 10X intervals from 11+(>10X) to 121+.

-Then to summarise the coverage stats per cohort

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
RESULTSDIR=${PROJECTDIR}/results/qc_plots/depth
DPSCRIPTDIR=${PROJECTDIR}/scripts/coverage_qc
#Togenerate the cov_stats summary for the cohort
cd ${RESULTSDIR}
#Get the header from the first file
for f in `ls -1 *.covstats.tsv | head -n1`; do head -n1 ${f} | sed 's/^11/Sample\t11/' >cov_stats_summary.tsv; done
grep -v Mean *stats.tsv | sed 's/.covstats.tsv:/\t/' >> cov_stats_summary.tsv

# Generate the tables and plots for the coverage stats summary
cd ${PROJECTDIR}
Rscript ${DPSCRIPTDIR}/02_plot_cov_stats_summary_sortByMeanCov_mod.R
```
