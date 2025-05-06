# 2744_IVO_Cherry_angioma_WES
# A post-zygotic disruptive germline _NFATC1_ variant in a patient with segmental cherry angiomas

[![DOI]()]()

This repository serves as the central landing page for mulitple sub-analyses related to the manuscript:

> **_A post-zygotic disruptive germline NFATC1 variant in a patient with segmental cherry angiomas_**

## Overview

This repository contains a series of scripts that were used in generating the results and figures for the manuscript. The analysis is divided into several sections, each corresponding to a specific aspect of the study. The main analyses include:

- **Variant calling and annotation**: This section includes the scripts used for variant calling from whole exome sequencing data. It includes the generation of VCF files and the filtering of variants based on various criteria.
- **Assessment of _NFATC1_ expression on endothelial cells**: This section looks at the expression of _NFATC1_ in endothelial cells across organs, annotated cells on the skin dataset and  differential expression analysis between Endothelial cells on the skin by _NFATC1_ expression from the _Tabula Sapiens_ single cell datasets.  It includes the analysis of RNA-seq data and the generation of relevant plots.


## Dependencies

This analysis was conducted within a development container with docker and R (4.3.3) installed (`.devcontainer/devcontainer.json`). All scripts can therefore be run in a Github codespace or Vscode session with docker available. All R dependencies for this project are detailed within the project `renv.lock` file and can be installed by running `renv::restore()` in a R terminal.


## Contact 
- Martin Del Castillo Velasco-Herrera - (<mdc1@sanger.ac.uk>)
- Alexis Germán Murillo Carrasco - (<alexis.carrasco@hc.fm.usp.br>)
- David J Adams - (<da1@sanger.ac.uk>)