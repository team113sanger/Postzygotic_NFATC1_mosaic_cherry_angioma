# dbsnp data preparation for Annotation

The dbSNP data is large and must be preprocessed for its use with the MAF and QC scripts used. The steps to generate the SNP files used for somatic calling are described below.

## Required software
- `wget` to download the dbSNP files
- `bgzip` and `tabix` to compress and index the VCF files
- `bcftools` to manipulate the VCF files

## Download the dbSNP files

dbSNP v155 was used. We need the VCF and the assembly report to convert contig names to chromosomes. Assign a `PROJECTDIR` variable to the project directory where repository was cloned.

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
DBSNPDIR=${PROJECTDIR}/resources/dbsnp
cd ${DBSNPDIR}

wget https://ftp.ncbi.nlm.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt

# The VCF is not compressed with bgzip; must unzip and use bgzip and then use tabix to index
gunzip GCF_000001405.39.gz
bgzip gunzip GCF_000001405.39
tabix -p vcf GCF_000001405.39.gz
```

## Covert chromosome names and create a smaller SNP file

Get the chromosome associations
```bash
grep Sequence-Name GCF_000001405.39_GRCh38.p13_assembly_report.txt -A 10000| awk 'BEGIN{FS="\t"}{ print $7,$10 }' > GCF_000001405.39_GRCh38.p13_assembly_report.chrnames
```

Use bcftools annotate to rename the chromosomes and then remove some non-primary chromosomes that were changed to 'na'

```bash
bsub -q normal  -M2000 -R"select[mem>2000] rusage[mem=2000]" \
  -e rename.155.e -o rename.155.o -n 8 "module load bcftools/1.9; \
  bcftools annotate --rename-chrs GCF_000001405.39_GRCh38.p13_assembly_report.chrnames \
  --threads 8  GCF_000001405.39.gz | grep -v ^na | grep -v \
  "contig=<ID=na>\" > \
  dbSNP155.GRCh38.GCF_000001405.39.mod.vcf; \
  bgzip dbSNP155.GRCh38.GCF_000001405.39.mod.vcf; \
  tabix -p vcf dbSNP155.GRCh38.GCF_000001405.39.mod.vcf.gz"
```

## Common SNPs

Get common SNPs - SNPs with AF > 0.01 in 1000 Genomes are marked as common
```bash
echo -e "#CHROM\tPOS\tREF\tALT\tdbSNP" > dbSNP155_common.tsv
zcat dbSNP155.GRCh38.GCF_000001405.39.mod.vcf.gz |grep ";COMMON"|cut -f 1,2,3,4,5| awk '{print $0"\tCOMMON"}' >> dbSNP155_common.tsv

module load bgzip/1.18
bgzip dbSNP155_common.tsv; tabix --skip-lines 1 -b 2 -e 2 -s 1 dbSNP155_common.tsv.gz
```

# NOTE

The `dbSNP_common.tsv.gz` and `dbSNP155.GRCh38.GCF_000001405.39.mod.vcf.gz` files are too large




