#! /bin/bash
# THe script is aimed to be run as 
# run_cgpFlagCaVEManvv1.17.2.sh config_file_wFLAGS versionID tumour_PDID normal_PDID SAMPLE_Pair_name BAM_DIR VCF_Main_DIR
#bash ./run_cgpFlagCaVEManvv1.17.2.sh DERMATLAS_final_flag.vcf.config.ini v1 PD52540a PD52540b PD52540a_vs_D52540b /lustre/projects/BAMS/WES /lustre/projects/results/CAVEMAN

#INPUT parameters
config=${1}
version=${2}
tumour=${3}
normal=${4}
#Name of the sample pair but for the file i.e. no "-"  e.g. PD52540a_vs_D52540b
sample=${5}
BAMDIR=${6}
VCFDIR=${7}
#The name of the ouput folder
pair="${tumour}_vs_${normal}"
samplefolder="${tumour}"

# paths to the VCF files on lustre
export VCF_IN=${VCFDIR}/${samplefolder}/${sample}.muts.ids.mnv.vcf.gz
export VCF_OUT=${VCFDIR}/${samplefolder}/${sample}.smartphase.flag.vcf
# paths to the BAM/CRAM files on lustre
export MT_BAM=${BAMDIR}/${tumour}/${tumour}.sample.dupmarked.bam 
export WT_BAM=${BAMDIR}/${normal}/${normal}.sample.dupmarked.bam

####################################################################################
# Special environment variables required for cgpcavemanwrapper/1.18.2
####################################################################################
# set this to your lustre area
export WORK=${PWD}
export VCFDIR=${VCFDIR}
export REF_FILE=/lustre/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa.fai
export BED_LOC=/lustre/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/flagging/
export FLAG_FILE=/lustre/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/flag.to.vcf.convert.ini
export UNMATCHED_VCF_LOC=/lustre/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/unmatched_v5_vcf
export ANNOT_BED_LOC=/lustre/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/vagrent/e103

export SINGULARITY_BINDPATH="/lustre,/nfs/users/nfs_m/mdc1,${WORK},"

echo " module load cgpcavemanwrapper/1.17.2 ;  
cgpFlagCaVEMan.pl \
    --input ${VCF_IN} \
    --outFile ${VCF_OUT} \
    --species Human \
    --reference ${REF_FILE} \
    --studyType WXS \
    --tumBam ${MT_BAM} \
    --normBam ${WT_BAM} \
    --bedFileLoc ${BED_LOC} \
    --unmatchedVCFLoc ${UNMATCHED_VCF_LOC} \
    --annoBedLoc ${ANNOT_BED_LOC} \
    --flagConfig ${config} \
    --flagToVcfConfig ${FLAG_FILE}

"
module load cgpcavemanwrapper/1.17.2

cgpFlagCaVEMan.pl \
    --input ${VCF_IN} \
    --outFile ${VCF_OUT} \
    --species Human \
    --reference ${REF_FILE} \
    --studyType WXS \
    --tumBam ${MT_BAM} \
    --normBam ${WT_BAM} \
    --bedFileLoc ${BED_LOC} \
    --unmatchedVCFLoc ${UNMATCHED_VCF_LOC} \
    --annoBedLoc ${ANNOT_BED_LOC} \
    --flagConfig ${config} \
    --flagToVcfConfig ${FLAG_FILE}


