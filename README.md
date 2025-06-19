# 2744_IVO_Cherry_angioma_WES
# A post-zygotic disruptive germline _NFATC1_ variant in a patient with segmental cherry angiomas

[![DOI]()]()

This repository serves as the central landing page for mulitple sub-analyses related to the manuscript:

> **_A post-zygotic disruptive germline NFATC1 variant in a patient with segmental cherry angiomas_**

## Overview

This repository contains a series of scripts that were used in generating the results and figures for the manuscript. The analysis is divided into several sections, each corresponding to a specific aspect of the study. The main analyses include:

- **Variant calling and annotation**: This section includes the scripts used for variant calling from whole exome sequencing data. It includes the generation of VCF files and the filtering of variants based on various criteria. To see [README](documentation/Somatic_Variant_calling.md) file in the `documentation` folder for more details on how to reproduce the analysis.
- **Assessment of _NFATC1_ expression on endothelial cells**: This section looks at the expression of _NFATC1_ in endothelial cells across organs, annotated cells on the skin dataset and  differential expression analysis between Endothelial cells on the skin by _NFATC1_ expression from the _Tabula Sapiens_ single cell datasets.  It includes the analyses and the generation of relevant plots for the figures. To reproduce the analysis, please refer to the [README](documentation/Tabula_sapiens_NFATC1_exp_analysis.md) file in the `documentation` folder.


## Dependencies

Users can reproduce the scripts bundled in this package by creating a Github codespace using the .devcontainer.json file that is bundled with these scripts. 

This analysis was conducted within a development container with docker and R (4.3.3) installed (`.devcontainer/devcontainer.json`). All scripts can therefore be run in a Github codespace or Vscode session with docker available. All R dependencies for this project are detailed within the project `renv.lock` file and can be installed by running `renv::restore()` in a R terminal.

### Software dependencies

Analyses were performed using a combination of R, Perl and shell scripts.

- All R scripts were run using `R v4.3.3` and the package dependencies for each analysis are detailed within the `renv.lock`.  
- Perl scripts were run using Perl version `v5.38.0`
- The following software was used:
  - `samtools` version`v1.14` [**here**](https://github.com/samtools/samtools)
  - `bwa-mem` version `0.7.17` [**here**](https://github.com/lh3/bwa)
  - `XenofilteR` version `1.6` [**here**](https://github.com/NKI-GCF/XenofilteR/releases/tag/v1.6)
  - `bcftools` version `1.9` [**here**](https://github.com/samtools/bcftools/)
  - `tabix` version `1.9` [**here**](https://github.com/samtools/tabix/)
  - `CaVEMan` version `1.18.2` [**here**](https://github.com/cancerit/CaVEMan)
  - `cgpCaVEManwrapper` version `1.18.2` [**here**](https://github.com/cancerit/cgpCaVEManWrapper)
  - `Smart-Phase` version `0.1.8`[casm docker image here](https://github.com/cancerit/CASM-Smart-Phase/tree/main)
  - `cgpPindel` version `3.11.0` [**here**](https://github.com/cancerit/cgpPindel)
  - ENSEMBL VEP version `103`[**here**](http://feb2021.archive.ensembl.org/info/docs/tools/vep/index.html)
- The repositories with the scripts used for variant QC and VCF to MAF conversion can be found in the following links:
    - [**QC**](https://github.com/team113sanger/dermatlas_analysis_qc) `v0.4.3`
    - [**MAF**](https://github.com/team113sanger/dermatlas_analysis_maf) `v0.5.4` 


## Contact 
- Martin Del Castillo Velasco-Herrera - (<mdc1@sanger.ac.uk>)
- Alexis Germán Murillo Carrasco - (<alexis.carrasco@hc.fm.usp.br>)
- David J Adams - (<da1@sanger.ac.uk>)

