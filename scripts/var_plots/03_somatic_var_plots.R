#!/usr/bin/env Rscript
# A script to genereate plots for the somatic variants observed

#Get the current working directory
here::i_am("scripts/var_plots/03_somatic_var_plots.R")

library(here)
library(maftools)
library(logger)

# Load the MAF file
projdir <- here::here()
maf_dir<- file.path(projdir, "analysis", "variants_combined", "version1")
results_dir <- file.path(projdir, "results", "somatic_var_plots")
fig_dir <- file.path(projdir,"results", "figures")
#Check if the directories exist, if not create them
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings= FALSE)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
}

logger::log_info("Reading the MAF file with the somatic variants")
# Load the MAF file
chang_maf_url<-"https://figshare.com/ndownloader/files/55690916"
chang_maf_file<-file.path(projdir, "analysis", "variants_combined", "version1", "8117_2744-filtered_mutations_matched_allTum_keep.maf")
chang_clin<-file.path(projdir, "metadata", "8117_2744_metadata.tsv")

#readRDS(url(tbs_skin_data_url))
chang_maf<-read.maf(maf = chang_maf_file,
                    clinicalData = chang_clin)

# Set colours for tissue type 
# #C8403D #B9975BFF
pcol<-c( "#C8403D","#B9975BFF")
names(pcol)<- unique(chang_maf@clinical.data$Tissue_type)
pcol<-list(Tissue_type=pcol)

oncoplot_fname<-file.path(results_dir, "Cherry_angioma_HS_som_mut_oncoplot.pdf")
pdf(file =oncoplot_fname ,width = 8, height = 8)
oncoplot(chang_maf,
         top = 20,
         altered=F,
         showTitle = F,
         showTumorSampleBarcodes = T,
         drawColBar = FALSE,
         SampleNamefontSize = 1,
         additionalFeatureCex = 0.9,
         #sampleOrder = c(), # This is to ensure the order as the tissue section 
         annotationFontSize = 1.6,
         legendFontSize = 1.6, 
         barcode_mar = 5,
         gene_mar=6,
         clinicalFeatures = 'Tissue_type',
         sortByAnnotation = FALSE,
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = pcol, # This is to ensure the colours match between slide and the samples
         draw_titv = F)
dev.off()


loli_fname<-file.path(results_dir, "NFATC1_loliplot.pdf")
pdf(file =loli_fname ,width = 7, height = 4)
lollipopPlot(
  maf = chang_maf,
  gene = 'NFATC1',
  AACol = 'HGVSp',
  showMutationRate = FALSE,
  labelPos = 741
)
dev.off()