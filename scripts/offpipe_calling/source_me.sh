#!/bin/bash

module purge

module load cgpcavemanwrapper/1.17.2
module load cgppindel/3.5.0


module load bcftools-1.9/python-3.11.6
module load ensembl_vep/103.1

#Baitset file
HUM_PAD_MERGED_BAITSET=/lustre/resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed

