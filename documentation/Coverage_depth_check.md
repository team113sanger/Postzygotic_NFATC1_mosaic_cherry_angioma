# Coverage depth check for somatic variant calling	

## Overview

To analyse the evenness of coverage depth of the samples, we used a file that contained the coordinates of the baits used for exome capture.  liftover to GRCh38 followed by the use of the `samtools depth` command to calculate the coverage depth of each sample. The output was then processed to generate a summary of the coverage depth across all samples.

Samples with a coverage depth of >20X coverage across at least 80% of the baits used were considered for somatic variant calling. 

## Dependencies

## Required software

- `CrossMap` version `0.5.4>` [**here**](http://crossmap.sourceforge.net/)

### Required data

- The bait set used for the exome capture was the Agilent SureSelect Human All Exon V6+UTR, the file was obtained from Agilent and obtained  (GRCh38 liftover). File [`resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed`](resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed) 
- The chain file for the liftOver from hg19 to GRCh38 was downloaded from UCSC [`hg19ToHg38.over.chain.gz`](resources/baits/chainset/hg19ToHg38.over.chain.gz)

-Since I Received the bait set information from Agilent SureSelectV6+UTRs bait set (S07604624_hs_hg19.zip)  located here:
${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_hg19

-I downloaded the UCSC chain file to be able to do the liftOver from hg19 to GRCh38  using CrossMaP( http://crossmap.sourceforge.net/  )
- The Chain file was downloaded from (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/ ) and the file was [hg19ToHg38.over.chain.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)
```bash
cd ${PROJECTDIR}/resources/baits/chainset
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

-Then to submit the jobs for the liftover of the BED files I did the following:
```bash
cd  ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_GRCh38
mkdir ./logs
bsub -q normal -o ./logs/cmap1_%J.o -e ./logs/cmap1_%J.e -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -n 6 '/software/team113/users/mdc1/python/bin/CrossMap.py bed ${PROJECTDIR}/resources/baits/chainset/hg19ToHg38.over.chain.gz ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_hg19/S07604624_Covered.bed ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_GRCh38/S07604624_Covered.bed '

bsub -q normal -o ./logs/cmap1_%J.o -e ./logs/cmap1_%J.e -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -n 6 '/software/team113/users/mdc1/python/bin/CrossMap.py bed ${PROJECTDIR}/resources/baits/chainset/hg19ToHg38.over.chain.gz ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_hg19/S07604624_AllTracks.bed ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_GRCh38/S07604624_AllTracks.bed'

bsub -q normal -o ./logs/cmap1_%J.o -e ./logs/cmap1_%J.e -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -n 6 '/software/team113/users/mdc1/python/bin/CrossMap.py bed ${PROJECTDIR}/resources/baits/chainset/hg19ToHg38.over.chain.gz ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_hg19/S07604624_Padded.bed ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_GRCh38/S07604624_Padded.bed '

bsub -q normal -o ./logs/cmap1_%J.o -e ./logs/cmap1_%J.e -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -n 6 'CrossMap.py bed ${PROJECTDIR}/resources/baits/chainset/hg19ToHg38.over.chain.gz ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_hg19/S07604624_Regions.bed ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_GRCh38/S07604624_Regions.bed'

cp ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_hg19/S07604624_Targets.txt ${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_GRCh38/
```

-The summary of the files  NUMBER of sequences that liftover properly and the ones that are missing from the important files present is :

```
mdc1@farm5-head2:${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_GRCh38$ wc -l S07604624_Padded.bed
299098 S07604624_Padded.bed
(base) mdc1@farm5-head2:${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_GRCh38$ wc -l S07604624_Padded.bed.unmap
20 S07604624_Padded.bed.unmap
(base) mdc1@farm5-head2:${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_GRCh38$ wc -l S07604624_Covered.bed
297506 S07604624_Covered.bed
(base) mdc1@farm5-head2:${PROJECTDIR}/resources/baits/Agilent_SureSelectV6_plus_UTRs_GRCh38$ wc -l S07604624_Covered.bed.unmap
23 S07604624_Covered.bed.unmap
```
- Edu has confirmed in an email  that  the samples were extracted from fresh frozen tissues. 
-I have a started the  somatic calling of the tumour samples with Caveman and Pindel 
-I have issued a sent a ticket [Hinxton 739371] Asking if there is a normal sample that  can be used with the GRCh38 ENS 103 reference


## Coverage Depth analysis
-To do the Coverage Depth analysis I have run the following:
TO CALCULATE THE  samples DEPTH for the entire COHORT
-To make a copy Links to the bam files I ran the following 
```bash
PROJECTDIR=/lustre/scratch117/casm/team113/projects/2744_IVO_Cherry_angioma_WES
cp -L /nfs/cancer_ref01/nst_links/live/2744/*/*.sample.dupmarked.bam $PROJECTDIR/BAMS/
#Then to Copy the indexes I ran the following
PROJECTDIR=/lustre/scratch117/casm/team113/projects/2744_IVO_Cherry_angioma_WES
cp -L /nfs/cancer_ref01/nst_links/live/2744/*/*.sample.dupmarked.bam.bai $PROJECTDIR/BAMS/
```

-To get the  copy of the liftovered baitset required I did the following: (method described in confluence here:
https://confluence.sanger.ac.uk/display/CIH/retrieving+the+baitset+used+by+a+project )
```bash
module load labProjectAdmin/2.3.2
PROJECTDIR=/lustre/scratch117/casm/team113/projects/2744_IVO_Cherry_angioma_WES
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
-Then to  summarise the coverage stats per sample
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

module load R/4.1.0; Rscript scripts/coverage_qc/02_plot_cov_stats_summary_sortByMeanCov_mod.R
```
