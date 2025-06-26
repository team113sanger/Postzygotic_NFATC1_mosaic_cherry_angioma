# Coverage depth check for somatic variant calling	

## Overview

To analyse the evenness of coverage depth of the samples, we used a file that contained the coordinates of the baits used for exome capture. Coordinates from hg19 were lifovered to GRCh38. Following this, coverage observed was performed using `samtools depth` command to calculate the coverage depth of each sample. The output was then processed to generate a summary of the coverage depth across all samples.

Samples with a coverage depth of >20X coverage across at least 80% of the baits used were considered for somatic variant calling. 

## Dependencies

## Required software

- Perl scripts were run using Perl version `v5.38.0`[**here**](https://www.perl.org/)
- `CrossMap` version `0.5.4>` [**here**](http://crossmap.sourceforge.net/)
	- The Chain file was downloaded from (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/ ) and the file was [hg19ToHg38.over.chain.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)
- `samtools` version `1.14` [**here**](https://www.htslib.org/)
- `R` version `4.3.3` [**here**](https://www.r-project.org/)

### Required datasets

- The bait set used for the exome capture was the Agilent SureSelect Human All Exon V6+UTR, the file was obtained from Agilent (file S07604624_hs_hg19.zip) and obtained  (liftover from hg19 to GRCh38 coordinates was done). See file [`resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed`](resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed) 
- The chain file for the liftOver from hg19 to GRCh38 used with CrossMap was downloaded from UCSC [`hg19ToHg38.over.chain.gz`](resources/baits/chainset/hg19ToHg38.over.chain.gz)

## Coverage Depth analysis

The Coverage Depth analysis was run under internal cluster environment. It is required to replace the `PROJECTDIR` shell environment variable with the path of the   I have run the following:

- To make a copy Links to the bam files I ran the following 
```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
cp -L /nfs/cancer_ref01/nst_links/live/2744/*/*.sample.dupmarked.bam $PROJECTDIR/BAMS/
#Then to Copy the indexes I ran the following
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
cp -L /nfs/cancer_ref01/nst_links/live/2744/*/*.sample.dupmarked.bam.bai $PROJECTDIR/BAMS/
```

-To get the  copy of the liftovered baitset required I did the following: (method described in confluence here:
https://confluence.sanger.ac.uk/display/CIH/retrieving+the+baitset+used+by+a+project )
```bash
module load labProjectAdmin/2.3.2
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
PROJECTDIR=/lustre/scratch125/casm/teams/team113/users/mdc1
mkdir -p $PROJECTDIR/required_datasets
cd $PROJECTDIR/required_datasets
#This is the "SureSelect_Human_All_Exon_V6+UTR (GRCh38 liftover)" baitset ID in canapps 177
get_baitset.pl -i 177
mv SureSelect_Human_All_Exon_V6+UTR\ \(GRCh38\ liftover\)_177.bed SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed
```

-To make a job array to calculate the depth with Samtools 
```bash
mkdir -p $PROJECTDIR/results/qc_plots/depth
mkdir -p $PROJECTDIR/logs/depth
cd $PROJECTDIR/results/qc_plots/depth
#To create the sample list file 
ls -1 $PROJECTDIR/BAMS/*.bam >sample_bams.fofn
ls -1 $PROJECTDIR/BAMS/*.bam | cut -f 9 -d "/" | sed 's/\.sample\.dupmarked\.bam//' >sample.list
#To create the jobs for submission of sample depth calculation
for f in `cat sample.list`; do bsub -e $PROJECTDIR/logs/depth/depth.$f.e -o $PROJECTDIR/logs/depth/depth.$f.o -n 2 -M2000 -R"select[mem>2000] rusage[mem=2000]" "module load samtools/1.14; samtools depth -a -b $PROJECTDIR/required_datasets/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed -o $f.depth.tsv -J -s -@ 2 $PROJECTDIR/BAMS/$f.sample.dupmarked.bam"; done

```
- Then to  summarise the coverage stats per sample

```bash
mkdir -p $PROJECTDIR/scripts
cp /lustre/scratch119/casm/team113da/projects/6500_DERMATLAS_SG_basal_cell_adenoma_and_adenocarcinoma_WES/SCRIPTS/DEPTH/* $PROJECTDIR/scripts
cd $PROJECTDIR/qc_plots/depth
for f in `cat sample.list`; do bsub -e $PROJECTDIR/logs/depth/count.$f.e -o $PROJECTDIR/logs/depth/count.$f.o -M2000 -R"select[mem>2000] rusage[mem=2000]" -q small "$PROJECTDIR/scripts/count_region_coverage.pl $f.depth.tsv > $f.covstats.tsv"; done
```
-Then to summarise the coverage stats per cohort
```bash
#Togenerate the cov_stats summary for the cohort

#Get the header from the first file
for f in `ls -1 *.covstats.tsv | head -n1`; do head -n1 $f | sed 's/^11/Sample\t11/' >cov_stats_summary.tsv; done
grep -v Mean *stats.tsv | sed 's/.covstats.tsv:/\t/' >> cov_stats_summary.tsv

Rscript scripts/coverage_qc/02_plot_cov_stats_summary_sortByMeanCov_mod.R
```
