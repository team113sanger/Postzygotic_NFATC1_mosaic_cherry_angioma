#!/usr/bin/env Rscript
# A script to genereate plots for the somatic variants observed

#Get the current working directory
here::i_am("scripts/var_plots/03_somatic_var_plots.R")

library(here)
library(maftools)

# Load the MAF file
projdir <- here::here()
maf_dir
results_dir <- file.path(projdir, "results", "somatic_var_plots")
fig_dir <- file.path(projdir, "figures", "somatic_var_plots")
#Check if the directories exist, if not create them
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings= FALSE)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
}


laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 





