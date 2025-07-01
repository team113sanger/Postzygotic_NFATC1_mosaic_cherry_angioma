# Somatic variant calling parameters, annotation filtering and results

## Overview

This document describes the steps taken to identify somatic single nucleotide variants (SNVs) and short INDEL mutations present on the whole exome of the biopsied samples. This was done using CaVEMAN and Pindel, followed by variant consequence prediction with VEP.  All samples were checked for evenness of coverage before the variant calling was performed (see [documentation/Coverage_depth_check.md](./Coverage_depth_check.md)).   

The somatic calling, flagging and annotation was performed with an internal pipeline. The steps, and parameters used for the variant calling are described below including some example scripts. However, **Filtering, MAF conversion and plotting** were performed using the commands and scripts mentioned below on their respective section

## Results

MAF file containing the summary of the results of the variant calling can be found in Figshare Project [here](https://figshare.com/projects/A_post-zygotic_disruptive_germline_NFATC1_variant_in_a_patient_with_segmental_cherry_angiomas/254267). 

Raw VCF files containing the flagged and VEP annotated variant calls can be downloaded from EGA study Accession **EGAS00001008212**


## Sample Pairs used for somatic calling 

The list of sample pairs used for variant calling using the somatic callers can be found in the file [`metadata/8117-biosample_manifest-completed.tsv`](../metadata/8117-biosample_manifest-completed.tsv). 
 
## Required Environment variables and software

The following software is required to be installed and visible in the path before running the scripts:
- **R**: R `4.3.3`[**here**](https://cran.r-project.org/)
- **samtools**: `v1.14` [**here**](https://github.com/samtools/samtools)
- **bcftools**: `1.9` [**here**](https://github.com/samtools/bcftools/)
- **CaVEMan**: `1.15.1` [**here**](https://github.com/cancerit/CaVEMan)
- **cgpCaVEManwrapper**:`1.17.2` [**here**](https://github.com/cancerit/cgpCaVEManWrapper)
- **cgpPindel**: `v3.5.0` [**here**](https://github.com/cancerit/cgpPindel)
- ENSEMBL VEP version `103`[**here**](http://feb2021.archive.ensembl.org/info/docs/tools/vep/index.html)
- The repositories with the scripts used for variant QC and VCF to MAF conversion have been added as **submodules** to this repository:
    - [**QC**](https://github.com/team113sanger/dermatlas_analysis_qc) `v0.4.3`
    - [**MAF**](https://github.com/team113sanger/dermatlas_analysis_maf) `v0.5.4` 
    - To clone the submodules used if not done when cloning the repository, run the following command:
```bash 
git submodule update --init --recursive
```

:warning: **IMPORTANT NOTES** :warning:
- The scripts below are written to be run in our internal servers, which have Ubuntu 22.04.5 and submit jobs using `lsf bsub`. They are written to show an example of the commands used, this is because the calling was performed with an internal pipeline that uses CaVEMan and Pindel which utilise the same singularity images can be made from the repositories mentioned above.  Paths and environment variables would need to be adjusted to run in a different environment.
- `cgpCaVEManwrapper` and `cgpPindel` were downloaded and used as singularity images within our internal pipelines using `bpipe`[https://docs.bpipe.org/](https://docs.bpipe.org/). 
- bpipe running parameters used for CaVEMan and Pindel can be found in the [bpipe_CaVEMan_Pindel_running_params.md](./bpipe_CaVEMan_Pindel_running_params.md) file.
- Path, or scripts module call modifications may need to be made to run in a different environment.

## Somatic variant calling

The variant calling was performed using `CaVEMan` and `cgpPindel` somatic callers. To see the samples pairs used as tumour-normal for the calling see [metadata/8117-biosample_manifest-completed.tsv](../metadata/8117-biosample_manifest-completed.tsv). Metadata for the samples can be located in the[metadata/8117_2744_metadata.tsv](../metadata/8117_2744_metadata.tsv) table. Functional annotation was done with ENSEMBL v103 Variant Effect Predictor (VEP).  All jobs of these steps were performed inside an internal pipeline however the running for all the callers and variant effect prediction were the same as the ones used below. 

####  **STEP 1- CAVEMAN v1.15.1 SNV calling  **

Running parameters for the Caveman calling can be found inside the [run_Cavemanwrapper_1.17.2.sh](../scripts/offpipe_calling/run_Cavemanwrapper_1.17.2.sh) script.

**NOTE**
Set the `PROJECTDIR` variable to the path where the repository was cloned into and run the following commands in the terminal:

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365
PREFIX="${STUDY}_${PROJECT}"

cd ${PROJECTDIR}/scripts

# BAM Directory for the Xenofilitered WES data
BMDIR=${PROJECTDIR}/bams/WES_xfilt/NOD_PDXV1

mkdir -p ${PROJECTDIR}/analysis/CAVEMAN
cd ${PROJECTDIR}/analysis/CAVEMAN
mkdir -p logs;
num=0; 
for Tumour_PDID in `cut -f 1 ${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25.tsv`;  do 
echo $num; 
let num=num+1; 
Normal_PDID=`grep ${Tumour_PDID} ${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25.tsv |cut -f 2`; 
bsub -e ./logs/caveman.${Tumour_PDID}-vs-$Normal_PDID.e -o ./logs/caveman.${Tumour_PDID}-vs-$Normal_PDID.o \
-J"cavemanrun[$num]" -n 7 -M40000 -R"select[mem>40000] rusage[mem=40000] span[hosts=1]"\
 -q long "bash ${PROJECTDIR}/scripts/offpipe_calling/run_Cavemanwrapper_1.17.2.sh ${Tumour_PDID:?unset} ${Normal_PDID:?unset} 6 ${BMDIR:?unset} "; 
done

```

#### **STEP 2- Call/Identify MNVs using SmartPhase 0.1.8 from unflagged CaVEMan**

There are 3 steps to call MNVs from unflagged CaVEMan calls using smart-phase:
- generate-bed (to get the adjacent MNVs)
- smart-phase (to phase the MNVs)
- merge-mnvs (adjust the original VCF to include MNVs)

##### A) generate -the bed files. From an unflagged CaVEMan VCF, create a BED file containing adjacent SNVs to check with SmartPhase.

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

# BAMDIR
BMDIR=${PROJECTDIR}/bams/WES_xfilt/NOD_PDXV1
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN

cd ${PROJECTDIR}/analysis/CAVEMAN
mkdir -p logs;
num=0; 
for Tumour_PDID in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
Normal_PDID=`grep ${Tumour_PDID} ${SAMPLETSV} |cut -f 2`;
PAIR=${Tumour_PDID}"_vs_"$Normal_PDID; 
bsub -e ./logs/smartphase.genbed.${Tumour_PDID}-vs-$Normal_PDID.e -o ./logs/smartphase.genbed.${Tumour_PDID}-vs-$Normal_PDID.o \
-J"smartphase_gbed[$num]" -n 2 -M8000 -R"select[mem>8000] rusage[mem=8000] span[hosts=1]"\
 -q normal "bash ${PROJECTDIR}/scripts/offpipe_calling/run_mnv_casmsmartphase_generate_bed.sh ${Tumour_PDID:?unset} ${Normal_PDID:?unset} ${PAIR:?unset} ${CAVEMANDIR:?unset} "; 
done

```

##### B) Run smart-phase -to phase the MNV

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

# BAMDIR
BMDIR=${PROJECTDIR}/bams/WES_xfilt/NOD_PDXV1
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN

cd ${PROJECTDIR}/analysis/CAVEMAN
mkdir -p logs;
num=0; 
for Tumour_PDID in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
Normal_PDID=`grep ${Tumour_PDID} ${SAMPLETSV} |cut -f 2`;
PAIR=${Tumour_PDID}"_vs_"$Normal_PDID;
bsub -e ./logs/smartphase.phase.${Tumour_PDID}-vs-$Normal_PDID.e -o ./logs/smartphase.phase.${Tumour_PDID}-vs-$Normal_PDID.o \
-J"smartphase_gbed[$num]" -n 2 -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]"\
 -q normal "bash ${PROJECTDIR}/scripts/offpipe_calling/run_mnv_casmsmartphase0.1.8.sh ${Tumour_PDID:?unset} ${Normal_PDID:?unset} ${PAIR:?unset} ${CAVEMANDIR:?unset} ${BMDIR:?unset} "; 
done
```

##### C) merge-mnvs -to merge the MNVs in the VCF

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

# BAMDIR
BMDIR=${PROJECTDIR}/bams/WES_xfilt/NOD_PDXV1
# Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
#Directory where the CaVEMan results are
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN

cd ${PROJECTDIR}/analysis/CAVEMAN

mkdir -p logs;
num=0; 
for Tumour_PDID in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
Normal_PDID=`grep ${Tumour_PDID} ${SAMPLETSV} |cut -f 2`;
PAIR=${Tumour_PDID}"_vs_"$Normal_PDID;
bsub -e ./logs/smartphase.mergemnv.${Tumour_PDID}-vs-$Normal_PDID.e -o ./logs/smartphase.mergemnv.${Tumour_PDID}-vs-$Normal_PDID.o \
-J"smartphase_mergmnv[$num]" -n 2 -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]"\
 -q normal "bash ${PROJECTDIR}/scripts/offpipe_calling/run_mnv_casmsmartphase0.1.8_mergemnvs.sh ${Tumour_PDID:?unset} ${Normal_PDID:?unset} ${PAIR:?unset} ${CAVEMANDIR:?unset} "; 
done
```

#### STEP 3- CGP CAVEMAN v18.2 flagging - cgpFlagCaVEMan

For this step the flagging .ini files are required, these contain the parameters for the flagging of the variants. A copy of these are shared with this repository in the resources folder. The files are:
- [**flag.to.vcf.convert.ini**](../resources/caveman/flag.to.vcf.convert.ini) -  contains conversions of flags to FLAG ID
- [**flag.vcf.config.ini**](../resources/caveman/flag.vcf.config.ini) - specifies which flags are cutoffs to use for WXS, WGS, etc. 

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

# BAMDIR
BMDIR=${PROJECTDIR}/bams/WES_xfilt/NOD_PDXV1
# Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
# Caveman directory
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN
CAVECONFIG=/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/flag.vcf.config.ini

cd ${PROJECTDIR}/analysis/CAVEMAN
num=0; 
for Tumour_PDID in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
Normal_PDID=`grep ${Tumour_PDID} ${SAMPLETSV} |cut -f 2`;
PAIR=${Tumour_PDID}"_vs_"$Normal_PDID; 
bsub -e ./logs/caveman.cgpFlag.${Tumour_PDID}-vs-$Normal_PDID.e -o ./logs/caveman.cgpFlag.${Tumour_PDID}-vs-$Normal_PDID.o \
-J"cavemanflag[$num]" -n 4 -M8000 -R"select[mem>8000] rusage[mem=8000] span[hosts=1]"\
 -q long "bash ${PROJECTDIR}/scripts/offpipe_calling/run_cgpFlagCaVEManvv1.15.1_postprocessingmnv.sh ${CAVECONFIG:?unset} v1 ${Tumour_PDID} $Normal_PDID ${PAIR:?unset} ${BMDIR:?unset} ${CAVEMANDIR:?unset} "; 
done
```

Then the VCFs were compressed and indexed using the following commands:

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

# Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
# Caveman directory
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN
cd ${CAVEMANDIR}
cat ${SAMPLETSV:?unset} | cut -f 1,2 | sed 's/\t/_vs_/g' >samples.list

module load bcftools-1.9/python-3.11.6
find ${PROJECTDIR}/analysis/CAVEMAN/*/* -name "*_DNA.smartphase.flag.vcf" -type f -exec bgzip {} \;
find ${PROJECTDIR}/analysis/CAVEMAN/*/* -name "*_DNA.smartphase.flag.vcf.gz" -exec tabix -p vcf {} \;

```
#### STEP 4- CGP CAVEMAN v1.17.4  VEP103 annotation - run_vep_caveman_vdermEns103.target.custom.grch38.sh

To predict the effect of the identified variants, ENSEMBLs' VEP v103 was used with additional custom annotations regarding their presence in external datasets such as:
    - gnomAD [`**v3.1**`](https://gnomad.broadinstitute.org/news/2019-10-gnomad-v3-0/)
    - dbSNP [`**v155**`](https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz) Instructions on how the file was parsed can be found [`resources/dbsnp`](../resources/dbsnp)
    - COSMIC [`**v97**`](https://cancer.sanger.ac.uk/cosmic/download/cosmic)
    - ClinVar [`**20220115**`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022/clinvar_20220115.vcf.gz)

To run the VEP annotation the following commands were used:

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

#Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
#Caveman directory
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN
cd ${CAVEMANDIR}
# List of caveman vcf files
ls -1 ${CAVEMANDIR}/*/*_DNA.smartphase.flag.vcf.gz >samples.list
#Run the VEP annotation
bash ${PROJECTDIR}/scripts/VEP/run_vep_pindel_Ens103.dermatlas.grch38.sh ./samples.list
#Then index the VCFs
module load bcftools-1.9/python-3.11.6
find ./*/* -name "*_DNA.smartphase.flag.vcf.gz" -exec tabix -p vcf {} \;
```

#### STEP 4.1- Filter the VCFs for PASSed variants located +-100bp of the bait target regions

Variants were filtered to keep only those that passed flagging and were called are within +-100bp of bait set targeted regions. BED file used provided. 

The following commands were used to filter the VCFs:
```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

#Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
#Baitset BED file paded100 bp
BEDFILE=${PROJECTDIR}/resources/baitset/GRCh38_WES5_canonical_pad100.merged.bed

#Caveman results directory
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN
cd ${CAVEMANDIR}
# List of caveman-SNV-MNV-Flagged-Annotated vcf files
ls -1 ${CAVEMANDIR}/*/*_DNA.smartphase.flag.vep.vcf.gz >samples.list
bash ${PROJECTDIR}/scripts/filterCaVEManPassTargetsFromVCF.sh ./samples.list ${BEDFILE:?unset} 

#Index the vcfs
module load bcftools-1.9/python-3.11.6
find ./*/* -name "*.smartphase.flag.vep.vcf.gz" -exec tabix -p vcf {} \;
```

Edit VCF file naming from Tumour_vs_Normal to Tumour only output to allow easy conversion to MAF 

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365
#Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
#Caveman results directory
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN

cd ${CAVEMANDIR}
# Then sort the naming so it is consistent with the one obtained by canapps pipeline
num=0; 
for f in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
g=`grep $f ${SAMPLETSV} |cut -f 2`;
mv ${CAVEMANDIR:?unset}/${f}/${f}_vs_${g}.smartphase.flag.vep.vcf.gz ${CAVEMANDIR:?unset}/${f}/${f}.smartphase.vep.vcf.gz ; 
mv ${CAVEMANDIR:?unset}/${f}/${f}_vs_${g}.smartphase.flag.vep.vcf.gz.tbi ${CAVEMANDIR:?unset}/${f}/${f}.smartphase.vep.vcf.gz.tbi ;
# TO copy the files to the directory they belonged for QC plotting
cp -R ${CAVEMANDIR:?unset}/${f} ${PROJECTDIR}/analysis/caveman_files/
done
```
#### STEP 6- Pindel calling   -CGP Pindel v3.11 calling 

The list of FLAGS for cgpPindel v3.11 is required to run the calling. A copy of the rules file is shared with this repository in the resources folder. The file is:
- [**WXS_Rules.lst**](../resources/pindel/WXS_Rules.lst) -

To run the calling 
```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

# BAMDIR
BMDIR=${PROJECTDIR}/bams/WES_xfilt/NOD_PDXV1
#Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
PINDELDIR=${PROJECTDIR}/analysis/PINDEL3.11
# path to germline indels bedfile
export GERM_INDEL=${PROJECTDIR}/resources/caveman/example/empty_germline_bed_for_caveman.bed
#Pindel Rules file - copy in ${PROJECTDIR}/resources/pindel
export CGPINDRULES=/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/pindel/WXS_Rules.lst

mkdir -p $PINDELDIR
cd $PINDELDIR
mkdir -p logs;
num=0; 
for f in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
g=`grep $f ${SAMPLETSV} |cut -f 2`; 
bsub -e logs/pindel.$f-vs-$g.e -o logs/pindel.$f-vs-$g.o \
-J"pindelrun[$num]" -n 6 -M40000 -R"select[mem>40000] rusage[mem=40000] \
span[hosts=1]" -q long "bash ${PROJECTDIR}/scripts/offpipe_calling/run_cgpPindel_3.11.sh ${f:?unset} ${g:?unset} 6 ${BMDIR:?unset} ${CGPINDRULES:?unset} "; 
done

```
#### STEP 6- CGP Pindel v3.11 calling  -POST-PROCESSING

To get the PASS calls within the targeted exome regions with +-100bp, the following commands were used:

- First the files were renames to remove the Tumour_vs_Normal format in the VCF file names 
```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

# BAMDIR
BMDIR=${PROJECTDIR}/bams/WES_xfilt/NOD_PDXV1
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
PINDELDIR=${PROJECTDIR}/analysis/PINDEL3.11

#Re name the samples from Tumour_vs_Normal to just Tumour
num=0; 
for f in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
g=`grep $f ${SAMPLETSV} |cut -f 2`;
mv ${PINDELDIR:?unset}/${f}/${f}_vs_${g}.flagged.vcf.gz ${PINDELDIR:?unset}/${f}/${f}.pindel.flagged.vcf.gz ; 
mv ${PINDELDIR:?unset}/${f}/${f}_vs_${g}.flagged.vcf.gz.tbi ${PINDELDIR:?unset}/${f}/${f}.pindel.flagged.vcf.gz.tbi ;
done
```
-To run the filter of PASSing variants within the targeted exome regions with 100bp, the following commands were used:

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

# BAMDIR
BMDIR=${PROJECTDIR}/bams/WES_xfilt/NOD_PDXV1
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
PINDELDIR=${PROJECTDIR}/analysis/PINDEL3.11
# path to germline indels bedfile
export GERM_INDEL=${PROJECTDIR}/resources/caveman/example/empty_germline_bed_for_caveman.bed
#Pindel Rules file - copy in resources/pindel
export CGPINDRULES=/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/pindel/WXS_Rules.lst
#Baitset BED file paded100 bp
HUM_PAD_MERGED_BAITSET=${PROJECTDIR}/resources/baitset/GRCh38_WES5_canonical_pad100.merged.bed

cd ${PINDELDIR}

# GET variants that are within the targeted exome regions with 100bp
ls -1 ${PINDELDIR:?unset}/*/*.flagged.vcf.gz >flag_samples_files.list
bash ${PROJECTDIR}/scripts/offpipe_calling/filterPindelPassTargetsFromVCF.sh ./flag_samples_files.list ${HUM_PAD_MERGED_BAITSET:?unset}

```

#### CGP Pindel v3.9 calling  -  VEP annotation

To predict the effect of the identified variants, ENSEMBLs' VEP v103 was used with additional custom annotations regarding their presence in external datasets such as:
    - gnomAD [`**v3.1**`](https://gnomad.broadinstitute.org/news/2019-10-gnomad-v3-0/)
    - dbSNP [`**v155**`](https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz) Instructions on how the file was parsed can be found [`resources/dbsnp`](../resources/dbsnp)
    - COSMIC [`**v97**`](https://cancer.sanger.ac.uk/cosmic/download/cosmic)
    - ClinVar [`**20220115**`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022/clinvar_20220115.vcf.gz)

To run the VEP annotation the following commands were used:

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

# BAMDIR
BMDIR=${PROJECTDIR}/bams/WES_xfilt/NOD_PDXV1
SAMPLETSV=${PROJECTDIR}/metadata/7688_3365_samplepairs_OMM25_all.tsv
PINDELDIR=${PROJECTDIR}/analysis/PINDEL3.11
# path to germline indels bedfile
export GERM_INDEL=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES/resources/caveman/example/empty_germline_bed_for_caveman.bed
#Pindel Rules file - copy in resources/pindel
export CGPINDRULES=/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/pindel/WXS_Rules.lst

cd $PINDELDIR

#to create the file with the sample list for analysis (taking the sample pairs names
ls -1 ${PINDELDIR:?unset}/*/*.flagged.target.pass.vcf.gz >flag_targ_pass_samples_files.list

bash ${PROJECTDIR}/scripts/VEP/run_vep_pindel_Ens103.dermatlas.grch38.sh ./flag_targ_pass_samples_files.list

#Then sort the naming Tumour_vs_Normal to Tumour only output to allow easy conversion to MAF
num=0; 
for f in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
g=`grep $f ${SAMPLETSV} |cut -f 2`;
mv ${PINDELDIR:?unset}/${f}/${f}.flagged.vcf.gz ${PINDELDIR:?unset}/${f}/${f}.pindel.flagged.vcf.gz ; 
mv ${PINDELDIR:?unset}/${f}/${f}.flagged.vcf.gz.tbi ${PINDELDIR:?unset}/${f}/${f}.pindel.flagged.vcf.gz.tbi ;
mv ${PINDELDIR:?unset}/${f}/${f}.flagged.target.pass.vep.vcf.gz ${PINDELDIR:?unset}/${f}/${f}.pindel.vep.vcf.gz ;
mv ${PINDELDIR:?unset}/${f}/${f}.flagged.target.pass.vep.vcf.gz.tbi ${PINDELDIR:?unset}/${f}/${f}.pindel.vep.vcf.gz.tbi ;
cp -R ${PINDELDIR:?unset}/${f} ${PROJECTDIR}/analysis/pindel_files/
done
```

### Additional annotations from dbSNP155 common variants

#### Adding annotations for dbSNP155 common variants  - CavEMan
```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365
STUDY_ID=3365

cd ${PROJECTDIR}/analysis
# Baitset BED file paded100 bp
BEDFILE=${PROJECTDIR}/resources/baitset/GRCh38_WES5_canonical_pad100.merged.bed
source ${PROJECTDIR}/scripts/QC/source_me.sh

cd ${PROJECTDIR}/analysis/caveman_files

# First get the on-target PASS variants from the CaVEMan file copied from nst_links (*smartphase.vep.vcf.gz)
# The output is *smartphase.vep.filt.vcf.gz
 
cd ${PROJECTDIR}/analysis/caveman_files

for f in */*.smartphase.vep.vcf.gz; do bash ${PROJECTDIR}/scripts/QC/select_vcf_pass_ontarget.sh $f $BEDFILE; done

# Now add dbSNP common annotations to the filtered VCF from above
for f in `dir -1`; do echo $f; bash ${PROJECTDIR}/scripts/QC/add_commonSNPs2vcf.sh -p ${PROJECTDIR} -o ./$f -v $f/$f.smartphase.vep.filt.vcf.gz; done
 
# This creates *snpflagged.vcf.gz files.
```

#### Adding annotations for dbSNP155 common variants  - Pindel

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365
STUDY_ID=3365

cd ${PROJECTDIR}/analysis

# First get the on-target PASS variants from the CaVEMan file copied from nst_links (*pindel.vep.vcf.gz)
# The output is *pindel.vep.filt.vcf.gz

# First get the on-target PASS variants from the Pindel files copied from nst_links (*pindel.vep.vcf.gz)
# The output is *pindel.vep.filt.vcf.gz
BEDFILE=${PROJECTDIR}/resources/baitset/GRCh38_WES5_canonical_pad100.merged.bed

source ${PROJECTDIR}/scripts/QC/source_me.sh
 
cd ${PROJECTDIR}/analysis/pindel_files
for f in */*pindel.vep.vcf.gz; do bash ${PROJECTDIR}/scripts/QC/select_vcf_pass_ontarget.sh $f ${BEDFILE}; done
 
# Now add dbSNP common annotations to the filtered VCF from above
# This creates *snpflagged.vcf.gz files.
 
for f in `dir -1 `; do echo $f; bash ${PROJECTDIR}/scripts/QC/add_commonSNPs2vcf.sh -p ${PROJECTDIR} -o ./$f -v $f/$f.pindel.vep.filt.vcf.gz; done
```
##### CREATE the sample lists with **QC/make_samplelists_from_manifest.pl**:

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365
STUDY_ID=3365
PREFIX="${STUDY}_${PROJECT}"

# Create the testing directory
source ${PROJECTDIR}/scripts/QC/source_me.sh
mkdir -p ${PROJECTDIR}/metadata
cd ${PROJECTDIR}/metadata

# USe the script  to generate the sample liests but using the rejected samples 
${PROJECTDIR}/scripts/QC/make_samplelists_from_manifest.pl --infile ${PROJECTDIR}/metadata/${STUDY}-biosample-manifest-completed.tsv --prefix ${STUDY}_${PROJECT}
# Optional parameter when needed not in the test --reject rejected
```

## Combine CaVEMan and Pindel files

To make MAF files and plots, the files containing the SNV, MNV, INDELs  were combined create a 'variants_combined' directory

```bash
# Create a subdirectory "version2" 
 mkdir -p ${PROJECTDIR}/analysis/variants_combined/version2
```

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

i=2
mkdir -p ${PROJECTDIR}/analysis/variants_combined/version${i}
mkdir -p ${PROJECTDIR}/analysis/variants_combined/release_v${i}
cd ${PROJECTDIR}/analysis/variants_combined/version${i}

source ${PROJECTDIR}/scripts/QC/source_me.sh

# If you have multiple independent tumours per patient and multiple tumours from the same tumour (duplicate samples), you can use "all" samples pulled from nst_links
dir -1  ${PROJECTDIR}/analysis/caveman_files/*/*.snpflagged.vcf.gz ${PROJECTDIR}/analysis/pindel_files/*/*snpflagged.vcf.gz | grep -f ${PROJECTDIR}/metadata/${STUDY}_${PROJECT}-analysed_all_tum.txt > caveman_pindel_vcfs_all.list
 
# To include all independent tumours per patient (no samples from the same tumour) - "independent"
dir -1  ${PROJECTDIR}/analysis/caveman_files/*/*.snpflagged.vcf.gz ${PROJECTDIR}/analysis/pindel_files/*/*.snpflagged.vcf.gz | grep -f ${PROJECTDIR}/metadata/${STUDY}_${PROJECT}-independent_tumours_all_tum.txt > caveman_pindel_vcfs_independent.list
 
# If you have multiple tumours per patient, exclude duplicate samples - "one per patient"
dir -1  ${PROJECTDIR}/analysis/caveman_files/*/*.snpflagged.vcf.gz ${PROJECTDIR}/analysis/pindel_files/*/*.snpflagged.vcf.gz | grep -f ${PROJECTDIR}/metadata/${STUDY}_${PROJECT}-one_tumour_per_patient_all_tum.txt > caveman_pindel_vcfs_onePerPatient.list
 

```
### Get the MAF files

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365
i=2

cd ${PROJECTDIR}/analysis/variants_combined/version${i}
source ${PROJECTDIR}/scripts/QC/source_me.sh
mkdir logs
# caveman_pindel_vcfs_${f}.list  - is your list of CaVEMan and Pindel VCFs
# caveman_pindel_${f}.maf  - is the suffix for the output MAF files
# filter2  - indicates filtering of indels based on VAF and size/type should be run
 
cd ${PROJECTDIR}/analysis/variants_combined/version${i}
 
for f in onePerPatient independent all; do
    mkdir $f
    cd $f
    mkdir -p logs
    bsub -e logs/qc.e -o logs/qc.o -M1200 -R"select[mem>1200] rusage[mem=1200]" -q normal \
    "bash ${PROJECTDIR}/scripts/QC/somatic_variants_qc.sh -l ../caveman_pindel_vcfs_${f}.list -m caveman_pindel_${f}.maf -s ${PROJECTDIR}/scripts -b GRCh38 -a gnomAD_AF -f filter2 -t ${PROJECTDIR}/resources/ensembl/dermatlas_noncanonical_transcripts_ens103.tsv"
    cd ..
done
```
OUTPUT:
The script outputs a list of MAF files and plots however relevant file for the next stage of the analysis is:
- [`analysis/variants_combined/version2/all/keep_caveman_pindel_all.maf`](../analysis/variants_combined/version2/all/keep_caveman_pindel_all.maf)


## Somatic Variant Filtering and plotting

Files  Adding annotations for dbSNP155 common variants  - CavEMan
```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
#It requires the canapps project ID instead
PROJECT=2744
cd ${PROJECTDIR}/analysis

BEDFILE=${PROJECTDIR}/resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed
source ${PROJECTDIR}/scripts/QC/source_me.sh

cd ${PROJECTDIR}/analysis/caveman_files
# First get the on-target PASS variants from the CaVEMan file copied from nst_links (*smartphase.vep.vcf.gz)
# The output is *smartphase.vep.filt.vcf.gz

for f in */*.caveman_c.vep.vcf.gz; do nice -n 5 bash ${PROJECTDIR:?unset}/scripts/QC/select_vcf_pass_ontarget.sh ${f} ${BEDFILE:?unset}; done

# Now add dbSNP common annotations to the filtered VCF from above
for f in `dir -1 |grep PD`; do echo $f; bash ${PROJECTDIR}/scripts/QC/add_commonSNPs2vcf.sh -p $PROJECTDIR -o ./$f -v $f/$f.caveman_c.vep.filt.vcf.gz; done
 
# This creates *snpflagged.vcf.gz files.
```

- Adding annotations for dbSNP155 common variants  - Pindel

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
#It requires the canapps project ID instead
PROJECT=2744
cd ${PROJECTDIR}/analysis

# First get the on-target PASS variants from the Pindel files copied (*pindel.vep.vcf.gz)
# The output is *pindel.vep.filt.vcf.gz
BEDFILE=${PROJECTDIR}/resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed
source ${PROJECTDIR}/scripts/QC/source_me.sh
 
cd ${PROJECTDIR}/analysis/pindel_files
for f in */*pindel.vep.vcf.gz; do bash ${PROJECTDIR}/scripts/QC/select_vcf_pass_ontarget.sh ${f} ${BEDFILE}; done
 
# Now add dbSNP common annotations to the filtered VCF from above
# This creates *snpflagged.vcf.gz files.
 for f in `dir -1 |grep PD`; do echo $f; bash ${PROJECTDIR}/scripts/QC/add_commonSNPs2vcf.sh -p $PROJECTDIR -o ./$f -v $f/$f.pindel.vep.filt.vcf.gz; done
```
- Make MAF files
```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744
i=1
# Variables
VCFLIST=${PROJECTDIR}/analysis/variants_combined/version${i}/caveman_pindel_vcfs_all.list
BUILD="GRCh38"
SAMPLELIST=${PROJECTDIR}/metadata/8117_2744-analysed_all.tsv
TEMP_MAF=${PROJECTDIR}/analysis/variants_combined/version${i}/caveman_pindel_pass_keep.maf
FINAL_MAF=${PROJECTDIR}/analysis/variants_combined/version${i}/${STUDY}_${PROJECT}-filtered_mutations_matched_allTum_keep.maf
SCRIPTDIR=${PROJECTDIR}/scripts/MAF

cd ${PROJECTDIR}/analysis/variants_combined/version${i}

# Generate a file with all the paths of the VCF files
dir -1  ${PROJECTDIR}/analysis/caveman_files/PD*/*snpflagged.vcf.gz ${PROJECTDIR}/analysis/pindel_files/PD*/*snpflagged.vcf.gz | grep -f ${PROJECTDIR}/metadata/${STUDY}_${PROJECT}-analysed_all_tum.txt > caveman_pindel_vcfs_all.list

# Transform the filtered VCFs into Combined MAFs
# Convert VCFs to maf
${SCRIPTDIR}/reformat_vcf2maf.pl --vcflist ${VCFLIST} --build ${BUILD} --pass --af_col "gnomAD_AF" --sample_list ${SAMPLELIST} --dbsnp_filter > ${TEMP_MAF}

#Get the keep variants - all variants in coding and splice sites regions (including synonymous variants)
head -n1 ${TEMP_MAF} > ${FINAL_MAF}
awk 'BEGIN{FS="\t"}{if($28~/keep/){print}}' ${TEMP_MAF} >> ${FINAL_MAF}

# Transform to Excel 
Rscript ${SCRIPTDIR}/maf2xlsx.R ${FINAL_MAF}

```
- Copy the files to the Results directory and  plot
```bash 
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744
i=1
# Variables
VCFLIST=${PROJECTDIR}/analysis/variants_combined/version${i}/caveman_pindel_vcfs_all.list
BUILD="GRCh38"
SAMPLELIST=${PROJECTDIR}/metadata/8117_2744-analysed_all.tsv
TEMP_MAF=${PROJECTDIR}/analysis/variants_combined/version${i}/caveman_pindel_pass_keep.maf
FINAL_MAF=${PROJECTDIR}/analysis/variants_combined/version${i}/${STUDY}_${PROJECT}-filtered_mutations_matched_allTum_keep.maf
SCRIPTDIR=${PROJECTDIR}/scripts/MAF

RESDIR
cd ${PROJECTDIR}/analysis/variants_combined/version${i}

Hugo_Symbol	Sample_ID	Protein_Change
NFATC1	PD54368a	S741Afs*17
```


## Plot the Variants obtained from the MAF files


**NOTE** If none of the other steps were run and you want to plot the variants from the MAF files, follow the steps below specified in the [`README.md`](../analysis/variants_combined/version1/README.md) file inside the `analysis/variants_combined/version1` directory. To get the MAF files from FigShare. 



