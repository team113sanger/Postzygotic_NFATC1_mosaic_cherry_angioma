
As example below we put the running parameters for sample "PD54368a" using `bpipe` pipeline for somatic calling with CaVEMan and Pindel. The parameters are used to run the analysis in our internal servers and may need to be adjusted to run in a different environment.

## Bpipe Canpipe running parameters for CaVEMan

### Bpipe running parameters for calling with CaVEMan
```bash
ALGORITHM => CaVEMan
COMMANDLINE => cd /lustre/canpipe/live/data/analysis/2744_3043054/CaVEMan; bpipe -r cgp_caveman.run
CONFIG_INI_FLAG => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/flag.vcf.config.ini
CONFIG_INI_VCF => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/flag.to.vcf.convert.ini
DB_TYPE => live
ENSEMBL_VER => 103
EXCLUDE_CONTIGS_FILE => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/ignore_contigs_caveman.txt
FIELD_TYPE => 3
FIELD_VAL => 3043054
FILE_TYPE => bam
FLAGGING_BED_DIR => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/flagging_v2
FLOWNAME => ANALYSIS CaVEMan
GENDER => XY
GENE_BUILD => 103
ID_ANALYSIS_PROC => 7595544
ID_INT_PROJECT => 2744
ID_INT_PROJECT_NORMAL => 2744
IGNORE_FILE => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed
INI_PATH => /software/CASM/modules/installs/canPipe/2.18.5/perl/config/canPipe.cfg.ini
JOBID => 660580
LOG4PERL => /software/CASM/modules/installs/canPipe/2.18.5/perl/config/log4perl_resultloader.conf
NORMAL_BAM => /lustre/canpipe/live/data/stage/2744_3043055/3043055.bam
NORMAL_CN_DEF => 2
NORMAL_CONTAMINATION => 0.00
NORMAL_ID => 3043055
NORMAL_NAME => PD54368b
NORMAL_PROTOCOL => WXS
NORM_SEQ_TYPE => WXS
NORM_STAGE => /lustre/canpipe/live/data/stage/2744_3043055
OUTDIR => /lustre/canpipe/live/data/analysis/2744_3043054/CaVEMan
PRIOR_PROB_MUT => 0.000006
PRIOR_PROB_SNP => 0.0001
REFERENCE => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa
REF_ROOT => /lustre/canpipe/live/ref/
SAMPLE_NAME => PD54368a
SCRATCH_ROOT => /lustre/canpipe/live/
SEQ_TYPE => WXS
SIF_CGPCAVEMANWRAPPER => /software/CASM/singularity/cgpcavemanwrapper/cgpcavemanwrapper_1.17.2.sif
SPECIES => Human
SPECIES_ASSEMBLY => GRCh38_full_analysis_set_plus_decoy_hla
TUMOUR_BAM => /lustre/canpipe/live/data/stage/2744_3043054/3043054.bam
TUMOUR_CN_DEF => 5
TUMOUR_ID => 3043054
TUMOUR_NAME => PD54368a
TUMOUR_PROTOCOL => WXS
TUM_STAGE => /lustre/canpipe/live/data/stage/2744_3043054
UNM_NORMAL_VCF_DIR => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/unmatched_v5
```
### Bpipe running parameters for CaVEMan flagging
```bash
ALGORITHM => CaVEMan_flagging
ANNOT_BED_DIR => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/vagrent/e${ENSEMBL_VER}/
COMMANDLINE => cd /lustre/canpipe/live/data/analysis/2744_3043054/CaVEMan_flagging; bpipe -r cgp_caveman_flagging.run
CONFIG_INI_FLAG => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/flag.vcf.config.ini
CONFIG_INI_VCF => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/flag.to.vcf.convert.ini
DB_TYPE => live
ENSEMBL_VER => 103
FIELD_TYPE => 3
FIELD_VAL => 3043054
FILE_TYPE => bam
FLAGGING_BED_DIR => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/flagging_v2
FLOWNAME => ANALYSIS CaVEMan_flagging
GENDER => XY
GENE_BUILD => 103
ID_ANALYSIS_PROC => 7601207
ID_INT_PROJECT => 2744
ID_INT_PROJECT_NORMAL => 2744
INI_PATH => /software/CASM/modules/installs/canPipe/2.18.5/perl/config/canPipe.cfg.ini
JOBID => 662212
LOG4PERL => /software/CASM/modules/installs/canPipe/2.18.5/perl/config/log4perl_resultloader.conf
NORMAL_BAM => /lustre/canpipe/live/data/stage/2744_3043055/3043055.bam
NORMAL_ID => 3043055
NORMAL_NAME => PD54368b
NORMAL_PROTOCOL => WXS
NORM_SEQ_TYPE => WXS
NORM_STAGE => /lustre/canpipe/live/data/stage/2744_3043055
OUTDIR => /lustre/canpipe/live/data/analysis/2744_3043054/CaVEMan_flagging
REFERENCE => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa
REF_ROOT => /lustre/canpipe/live/ref/
SAMPLE_NAME => PD54368a
SCRATCH_ROOT => /lustre/canpipe/live/
SEQ_TYPE => WXS
SIF_CGPCAVEMANWRAPPER => /software/CASM/singularity/cgpcavemanwrapper/cgpcavemanwrapper_1.17.2.sif
SPECIES => Human
SPECIES_ASSEMBLY => GRCh38_full_analysis_set_plus_decoy_hla
TUMOUR_BAM => /lustre/canpipe/live/data/stage/2744_3043054/3043054.bam
TUMOUR_ID => 3043054
TUMOUR_NAME => PD54368a
TUMOUR_PROTOCOL => WXS
TUM_STAGE => /lustre/canpipe/live/data/stage/2744_3043054
UNM_NORMAL_VCF_DIR => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/unmatched_v5
```
### Bpipe running parameters for CaVEMan  VEP
```bash
ALGORITHM => VEP_CaVEMan
COMMANDLINE => cd /lustre/canpipe/live/data/analysis/2744_3043054/VEP_CaVEMan; bpipe -r cgp_vep_caveman.run
DB_TYPE => live
ENSEMBL_VER => 103
FIELD_TYPE => 3
FIELD_VAL => 3043054
FILE_TYPE => bam
FLOWNAME => ANALYSIS VEP_CaVEMan
GENDER => XY
GENE_BUILD => 103
ID_ANALYSIS_PROC => 11798900
ID_INT_PROJECT => 2744
ID_INT_PROJECT_NORMAL => NA
ID_REFSET_GROUP => 698a3102-8681-48f3-b3d6-fd405bb3aacf
INI_PATH => /software/CASM/modules/installs/canPipe/2.53.5/perl/config/canPipe.cfg.ini
JOBID => 753951
LOG4PERL => /software/CASM/modules/installs/canPipe/2.53.5/perl/config/log4perl_resultloader.conf
NORMAL_BAM => NA/NA.bam
NORMAL_ID => NA
NORMAL_NAME => NA
NORMAL_PROTOCOL => WXS
NORM_SEQ_TYPE => NA
OUTDIR => /lustre/canpipe/live/data/analysis/2744_3043054/VEP_CaVEMan
REF_ROOT => /lustre/canpipe/live/ref/
SAMPLE_NAME => PD54368a
SCRATCH_ROOT => /lustre/canpipe/live/
SEQ_TYPE => WXS
SIF_ENSEMBL_VEP => /software/CASM/singularity/ensembl-vep/ensembl-vep_103.sif
SIF_PCAP_CORE => /software/CASM/singularity/pcap-core/pcap-core_5.6.1.sif
SPECIES => Human
SPECIES_ASSEMBLY => GRCh38_full_analysis_set_plus_decoy_hla
TUMOUR_BAM => /lustre/canpipe/live/data/stage/2744_3043054/3043054.bam
TUMOUR_ID => 3043054
TUMOUR_NAME => PD54368a
TUMOUR_PROTOCOL => WXS
TUM_STAGE => /lustre/canpipe/live/data/stage/2744_3043054
refBase => reference/
```
## Bpipe running parameters for Pindel

### Bpipe running parameters for somatic calling with Pindel
```bash
ALGORITHM => Pindel
BAD_LOCI => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz
COMMANDLINE => cd /lustre/canpipe/live/data/analysis/2744_3043054/Pindel; bpipe -r cgp_pindel.run
DB_TYPE => live
ENSEMBL_VER => 103
EXCLUDE_CONTIGS => chrUn%,HLA%,%_alt,%_random,chrM,chrEBV
EXCLUDE_CONTIGS_FILE => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/pindel/ignore_contigs_pindel.txt
FIELD_TYPE => 3
FIELD_VAL => 3043054
FILE_TYPE => bam
FILTER => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/pindel/WXS_Rules.lst
FLOWNAME => ANALYSIS Pindel
GENDER => XY
GENES => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/vagrent/e${ENSEMBL_VER}/codingexon_regions.indel.bed.gz
GENE_BUILD => 103
HIGH_DEPTH_REGIONS => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz
ID_ANALYSIS_PROC => 7595541
ID_INT_PROJECT => 2744
ID_INT_PROJECT_NORMAL => 2744
INI_PATH => /software/CASM/modules/installs/canPipe/2.18.5/perl/config/canPipe.cfg.ini
JOBID => 660576
LOG4PERL => /software/CASM/modules/installs/canPipe/2.18.5/perl/config/log4perl_resultloader.conf
NORMAL_BAM => /lustre/canpipe/live/data/stage/2744_3043055/3043055.bam
NORMAL_ID => 3043055
NORMAL_NAME => PD54368b
NORMAL_PROTOCOL => WXS
NORM_SEQ_TYPE => WXS
NORM_STAGE => /lustre/canpipe/live/data/stage/2744_3043055
OUTDIR => /lustre/canpipe/live/data/analysis/2744_3043054/Pindel
REFERENCE => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa
REF_ROOT => /lustre/canpipe/live/ref/
SAMPLE_NAME => PD54368a
SCRATCH_ROOT => /lustre/canpipe/live/
SEQ_TYPE => WXS
SIF_CGPPINDEL => /software/CASM/singularity/cgppindel/cgppindel_v3.5.0.sif
SIF_PCAP_CORE => /software/CASM/singularity/pcap-core/pcap-core_5.6.1.sif
SIF_VAFCORRECT => /software/CASM/singularity/vafcorrect/vafcorrect_5.7.0.sif
SIMREP => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/pindel/simpleRepeats.bed.gz
SOFT_FILTER => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/pindel/softRulesFragment.lst
SPECIES => Human
SPECIES_ASSEMBLY => GRCh38_full_analysis_set_plus_decoy_hla
TUMOUR_BAM => /lustre/canpipe/live/data/stage/2744_3043054/3043054.bam
TUMOUR_ID => 3043054
TUMOUR_NAME => PD54368a
TUMOUR_PROTOCOL => WXS
TUM_STAGE => /lustre/canpipe/live/data/stage/2744_3043054
UNM_NORMAL_VCF => /lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/pindel/pindel_np.v5.gff3.gz
```
### Bpipe running parameters for Pindel VEP annotation
```bash
ALGORITHM => VEP_Pindel
COMMANDLINE => cd /lustre/canpipe/live/data/analysis/2744_3043054/VEP_Pindel; bpipe -r cgp_vep_pindel.run
DB_TYPE => live
ENSEMBL_VER => 103
FIELD_TYPE => 3
FIELD_VAL => 3043054
FILE_TYPE => bam
FLOWNAME => ANALYSIS VEP_Pindel
GENDER => XY
GENE_BUILD => 103
ID_ANALYSIS_PROC => 11798898
ID_INT_PROJECT => 2744
ID_INT_PROJECT_NORMAL => NA
ID_REFSET_GROUP => 698a3102-8681-48f3-b3d6-fd405bb3aacf
INI_PATH => /software/CASM/modules/installs/canPipe/2.53.5/perl/config/canPipe.cfg.ini
JOBID => 753949
LOG4PERL => /software/CASM/modules/installs/canPipe/2.53.5/perl/config/log4perl_resultloader.conf
NORMAL_BAM => NA/NA.bam
NORMAL_ID => NA
NORMAL_NAME => NA
NORMAL_PROTOCOL => WXS
NORM_SEQ_TYPE => NA
OUTDIR => /lustre/canpipe/live/data/analysis/2744_3043054/VEP_Pindel
REF_ROOT => /lustre/canpipe/live/ref/
SAMPLE_NAME => PD54368a
SCRATCH_ROOT => /lustre/canpipe/live/
SEQ_TYPE => WXS
SIF_ENSEMBL_VEP => /software/CASM/singularity/ensembl-vep/ensembl-vep_103.sif
SIF_PCAP_CORE => /software/CASM/singularity/pcap-core/pcap-core_5.6.1.sif
SPECIES => Human
SPECIES_ASSEMBLY => GRCh38_full_analysis_set_plus_decoy_hla
TUMOUR_BAM => /lustre/canpipe/live/data/stage/2744_3043054/3043054.bam
TUMOUR_ID => 3043054
TUMOUR_NAME => PD54368a
TUMOUR_PROTOCOL => WXS
TUM_STAGE => /lustre/canpipe/live/data/stage/2744_3043054
refBase => reference/
```
