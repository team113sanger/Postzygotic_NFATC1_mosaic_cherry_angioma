#!/usr/bin/env Rscript
#Get the current working directory
here::i_am("scripts/tsapiens_analysis/04_tabula_sapiens_endothelial_cells_analysis.R")
# Tabula Sapiens - Skin analysis

library(here)
library(SeuratDisk)
library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(EnsDb.Hsapiens.v79)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(paletteer)
library(logger)

##################### Reading the required files  #######################
##Set work directory
projdir <- here::here()
tabdatadir<- file.path(projdir, "data", "tabula_sapiens")
setwd(projdir)
results_dir <- file.path(projdir, "results", "tsapiens_analysis")
fig_dir <- file.path(projdir, "results", "figures")
#Check if the results directory exists, if not create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
}
if(!dir.exists(fig_dir)){
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
}

## Get the data from the Tabula Sapiens - Skin data from CellxGene 
# Tabula Sapiens 1.0 Skin  URL https://datasets.cellxgene.cziscience.com/9405995f-b1c6-485c-8cc9-b80b9e870057.h5ad
# Download  .RDS with the Skin data from Tabula Sapiens 1.0 - ID 9405995f-b1c6-485c-8cc9-b80b9e870057.rds
# hd5ad files generated for TBS1.0 have issues in being converted to h5Seurat format to be loaded in Seurat
# Because of this, we are using the .rds file directly from the CellXFene repository
tbs_skin_data_url<-"https://datasets.cellxgene.cziscience.com/9405995f-b1c6-485c-8cc9-b80b9e870057.rds"
logger::log_info("Reading the Tabula Sapiens 1.0 skin file from the CellXgene repository")
logger::log_info(paste0("Reading : ", tbs_skin_data_url))
data_SKin <- readRDS(url(tbs_skin_data_url))

##Save first plot, dimplot
# With Colours more friendly to colourblind 
#FBE183 #F7D144 #F4C00D #FAAA05 #F78C09 #E35D2A #C8403D mast cell
 #A7373F #B14054 #D55474 #E26781 #E77988 #E78E97 #E3A0A5
 #C38AA3 #A8759F #A26296 #924E8A #723B79 #4F447E #2A6494
 #237E91 #299783 #5AAB6B #92C051
ident_colors <- colorRampPalette(paletteer::paletteer_d("MetBrewer::Signac",14 ))(length(unique(data_SKin$cell_ontology_class))) #Get colors for each cell type
ident_colors[11]<-"#C8403D" # Reassign Endothelial cell color
ident_colors[7]<-"#E26781" # Reassign Mast cell color
names(ident_colors) <- unique(data_SKin$cell_ontology_class) #Set names for each color
png(file.path(fig_dir, "Fig1C_TBS_skin_DimPlot_colblind.png"),width=7500,height=3000,res=600)
DimPlot(data_SKin, reduction = "umap",group.by = "cell_ontology_class", cols=ident_colors)+
#Reduce the size of the legend
  theme(legend.text=element_text(size=11))+
  labs(title = "Tabula Sapiens (Skin)")
dev.off()
pdf(file.path(fig_dir, "Fig1C_TBS_skin_DimPlot_colblind.pdf"),width=15,height=8)
DimPlot(data_SKin, reduction = "umap",group.by = "cell_ontology_class", cols=ident_colors)+
#Reduce the size of the legend
  theme(legend.text=element_text(size=11))+
  labs(title = "Tabula Sapiens (Skin)")
dev.off()



# Grey everything Else Endothelial cell only
gidentcol<-rep("lightgrey",length(unique(data_SKin$cell_ontology_class)))
gidentcol[11]<-"#C8403D" # Re-assign Endothelial cell color
names(gidentcol) <- unique(data_SKin$cell_ontology_class) #Set names for each color
png(file.path(results_dir, "TBS_skin_DimPlot_gendothelial.png"),width=7500,height=3000,res=600)
DimPlot(data_SKin, reduction = "umap",group.by = "cell_ontology_class", cols=gidentcol)+
#Reduce the size of the legend
  theme(legend.text=element_text(size=11))+
  labs(title = "Tabula Sapiens (Skin)")
dev.off()

##Set identity
data_SKin <- SetIdent(data_SKin, value = "cell_ontology_class")

##Sorting data by NFATC1 expression
a<-DotPlot(object = data_SKin, features = "ENSG00000131196")
a$data%>%arrange(-avg.exp)%>%select(id)%>%as.matrix()

data_SKin@active.ident <- factor(data_SKin@active.ident, 
                            levels=a$data%>%arrange(pct.exp)%>%select(id)%>%as.matrix()
)

logger::log_info("Plotting the expression of NFATC1 across all cell types")
png(file.path(fig_dir, "Fig1D_NFATC1_Exp_TS_skin_cells.png"),width=4500,height=3500,res=600)
DotPlot(object = data_SKin, features = "ENSG00000131196")+
  scale_x_discrete(labels=c( "NFATC1"))+xlab("")+
  theme(axis.text.x = element_text(face='italic'))
dev.off()
pdf(file.path(fig_dir, "Fig1D_NFATC1_Exp_TS_skin_cells.pdf"),width=8,height=7)
DotPlot(object = data_SKin, features = "ENSG00000131196")+
  scale_x_discrete(labels=c( "NFATC1"))+xlab("")+
  theme(axis.text.x = element_text(face='italic', size=15), 
        axis.text.y = element_text(size=13))
dev.off()


#plot with NFATC1 and VEGFR2 VEGFR1
logger::log_info("Plotting the expression of NFATC1, VEGFR1, VEGFR2 and MAPK1 (ERK) across all cell types")
png(file.path(results_dir, "TB_skin_NFATC1_VGFR1_2_MAPK1_Exp_TS_skin_cells.png"),width=4500,height=3500,res=600)
DotPlot(object = data_SKin, features = c("ENSG00000131196", "ENSG00000102755", "ENSG00000128052", "ENSG00000100030" ))+
  scale_x_discrete(labels=c( "NFATC1", "VEGFR1", "VEGFR2", "MAPK1"))+xlab("")+
  theme(axis.text.x = element_text(face='italic', angle=45, hjust=1))
dev.off()

########### Tabula Sapiens 1.0 - Endothelial cells All tissues
### Generate the plot of NFATC1 expression on endothelial cells across all tissues on the Tabula Sapiens v1.0
all_tissues_tbs_url<-"https://datasets.cellxgene.cziscience.com/4676160f-b0f6-4dea-a104-c2ef86c74674.rds"
logger::log_info("Reading the Tabula Sapiens 1.0 all organs RDS from the CellXgene repository")
logger::log_info(paste0("Reading : ", all_tissues_tbs_url))
tbs_v1_organs<-readRDS(url(all_tissues_tbs_url))

##Set identity of cells based on the tissue they belong to on the publication
tbs_v1_organs <- SetIdent(tbs_v1_organs, value = "tissue_in_publication")

##Sorting data by NFATC1
a<-DotPlot(object = tbs_v1_organs, features = "ENSG00000131196")
a$data%>%arrange(-avg.exp)%>%select(id)%>%as.matrix()
tbs_v1_organs@active.ident <- factor(tbs_v1_organs@active.ident, 
                            levels=a$data%>%arrange(pct.exp)%>%select(id)%>%as.matrix()
)

# Plotting the expression of NFATC1 across all tissues
png(file.path(fig_dir, "Fig1E_Exp_NFATC1_tissues_tbs1.png"),width=3000,height = 3500,res=600)
DotPlot(object = tbs_v1_organs, features = "ENSG00000131196")+
  scale_x_discrete(labels=c( "NFATC1"))+xlab("")+
  theme(axis.text.x = element_text(face='italic'))
dev.off()
pdf(file.path(fig_dir, "Fig1E_Exp_NFATC1_tissues_tbs1.pdf"),width=5,height = 7)
DotPlot(object = tbs_v1_organs, features = "ENSG00000131196")+
  scale_x_discrete(labels=c( "NFATC1"))+xlab("")+
  theme(axis.text.x = element_text(face='italic', size=15))
dev.off()


# Plotting the expression of NFATC1 across all tissues with VEGFR1, VEGFR2 and MAPK1 (ERK)
png(file.path(results_dir, "TS_WB_endo_Exp_NFATC1_VEGFR1_2_MAPK1_tissues_tbs1.png"),width=3000,height = 3500,res=600)
DotPlot(object = tbs_v1_organs, features = c("ENSG00000131196", "ENSG00000102755", "ENSG00000128052", "ENSG00000100030" ))+
  scale_x_discrete(labels=c( "NFATC1", "VEGFR1", "VEGFR2", "MAPK1"))+xlab("")+
  theme(axis.text.x = element_text(face='italic', angle=45, hjust=1))
dev.off()
rm(tbs_v1_organs)

##################### Tabula Sapiens 1.0 - Skin analysis #######################
#### Selecting endothelial cells
rownames(data_SKin@meta.data)[which(data_SKin@meta.data$cell_ontology_class=="endothelial cell")]
counts <- GetAssayData(data_SKin, slot="counts", assay="RNA") 
counts.sub <- counts[,which(data_SKin@meta.data$cell_ontology_class=="endothelial cell")] #only endothelial cells
metas.sub<-data_SKin@meta.data[which(data_SKin@meta.data$cell_ontology_class=="endothelial cell"),]
data_SKin2 <- CreateSeuratObject(counts=counts.sub,meta.data = metas.sub)

###Dichotomising dataset based on NFATC1 expression
logger::log_info("Dichotomising Tabula Sapiens Skin dataset based on NFATC1 expression")
med_NFATC1<- median(data_SKin2@assays$RNA$counts["ENSG00000131196",])

### Boxplot of log10 NFATC1 levels
png(file.path(results_dir, "Boxplot_NFATC1.png"),width=4500,height=3500,res=600)
boxplot(log10(data_SKin2@assays$RNA$counts["ENSG00000131196",]+1))
dev.off()

##Estimating cells expressing (or not) NFATC1
exp<-data.frame(NFATC1=data_SKin2@assays$RNA$counts["ENSG00000131196",])
exp<-exp%>%
  mutate(Group=ifelse(NFATC1==0,"Skin Endothelial cells NFATC1neg","Skin Endothelial cells NFATC1pos"))
table(exp$Group)
# expressing NFATC1 No expressing NFATC1 
# 510                  799 
exp$Group->data_SKin2@meta.data$groups

#Subset endothelial cells only
data_Skin3<-subset(data_SKin, cells=rownames(data_SKin@meta.data)[which(data_SKin@meta.data$cell_ontology_class=="endothelial cell")])
data_Skin3@meta.data$groups<-exp$Group
data_Skin3<- SetIdent(data_Skin3, value = "groups")

logger::log_info("Plotting the expression of NFATC1, VEGFR1, VEGFR2 and MAPK1 (ERK) by NFATC1 expression group")
#Plot  the expression NFATC1, VEGFR1, VEGFR2 and MAPK1 (ERK) FeaturePLot with NFATC1 and VEGFR2 VEGFR1
# Relevant for proliferation
# https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2020.599281/full
png(file.path(fig_dir, "Fig1H_TBS_skin_Endothelial_cells_by_NFATC_expgroup.png"),width=4500,height=3500,res=600)
DotPlot(object = data_Skin3, features = c("ENSG00000131196", "ENSG00000102755", "ENSG00000128052", "ENSG00000100030" ))+
  scale_x_discrete(labels=c( "NFATC1", "VEGFR1", "VEGFR2", "MAPK1"))+xlab("")+
  theme(axis.text.x = element_text(face='italic', angle=45, hjust=1))
dev.off()
pdf(file.path(fig_dir, "Fig1H_TBS_skin_Endothelial_cells_by_NFATC_expgroup.pdf"),width=7,height=6)
DotPlot(object = data_Skin3, features = c("ENSG00000131196", "ENSG00000102755", "ENSG00000128052", "ENSG00000100030" ))+
  scale_x_discrete(labels=c( "NFATC1", "VEGFR1", "VEGFR2", "MAPK1"))+xlab("")+
  theme(axis.text.x = element_text(face='italic', angle=45, hjust=1))
dev.off()

png(file.path(results_dir, "TBS_skin_Endocells_DimPlot_by_NFATC1_expgroup.png"),width=5500,height=3000,res=600)
DimPlot(data_Skin3, reduction = "umap", group.by = "groups")+
  labs(title = "Tabula Sapiens (Skin)")
dev.off()

###Calculating DESEQ manually
rm(list=c("a","exp","counts.sub","metas.sub","counts"))

####OTHER features
###Converting to a DESeq2 object
counts_ls <- list()
counts_ls <- data_SKin2@assays$RNA$counts+1

str(counts_ls)
# Plot the distribution of counts
png(file.path(results_dir, "Boxplot_counts_ls.png"),width=4500,height=3500,res=600)
boxplot(log10(rowSums(counts_ls)))
dev.off()

log10(sum(counts_ls["ENSG00000131196",])) #4.403 NFATC1 sum level 

## In this case, no exclusion of genes will be applied
counts_ls[which(log10(rowSums(counts_ls))>3.15),]->counts_ls ## 3.11 is log10(1309) the number of participant cells

logger::log_info("Performing differential expression analysis between NFATC1 expression groups using DESeq2")
dds <- DESeqDataSetFromMatrix(counts_ls, 
                              colData=data_SKin2@meta.data,
                              design = ~groups)

### running differential expression analysis
dds <- DESeq(dds)
res<-results(dds, tidy=TRUE)

############ Addition of the gene names to the results
dados<-res

##Converting ENSEMBL gene IDs to gene symbols and entrez IDs
logger::log_info("Genetting gene symbols and entrez IDs")
ensembl.genes <- dados$row
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID","ENTREZID"))
colnames(geneIDs1)<-c("Symbol","row","entrez")
dados<-merge(dados,geneIDs1,by="row")

logger::log_info("writing and plotting the differential expression results")
#Writin the differentially expressed genes  
write.csv(dados,file.path(results_dir,"DE_genes_by_NFATC1exp_results_Tidy.csv"), row.names=FALSE, quote=FALSE) ##Save the differential expression result


# Write a tbale with the DE results for NFATC1, VEGFR1, VEGFR2 and MAPK1 (ERK)
write.table(dados[dados$row %in% c("ENSG00000131196", "ENSG00000102755", "ENSG00000128052", "ENSG00000100030" ),], 
          file.path(results_dir, "DE_genes_by_NFATC1exp_results_Tidy_NFATC1_VEGFR1_2_MAPK1.tsv"), 
          row.names=FALSE, quote=FALSE, sep="\t") ##Save the differential expression result

############ Volcano plot
#dados$log2FoldChange*-1 -> dados$log2FoldChange
png(file.path(fig_dir, "Fig1F_VolcanoPlot_Tabula_endothelial.png"),width=4500,height = 3500,res=600)
EnhancedVolcano(dados,lab=dados$Symbol,y="padj",x="log2FoldChange",
                pCutoff = 1e-50,
                drawConnectors = T,
                title = substitute(paste(italic('NFATC1+ vs. NFATC1-'), "(subset: endothelial cells)")),
                subtitle = "p cutoff<1e-50 | FC cutoff > 2",
                gridlines.minor = FALSE,
                gridlines.major = FALSE,
                legendPosition = "bottom")
dev.off()
pdf(file.path(fig_dir, "Fig1F_VolcanoPlot_Tabula_endothelial.pdf"),width=10,height = 7)
EnhancedVolcano(dados,lab=dados$Symbol,y="padj",x="log2FoldChange",
                pCutoff = 1e-50,
                drawConnectors = T,
                title = substitute(paste(italic('NFATC1+ vs. NFATC1-'), "(subset: endothelial cells)")),
                subtitle = "p cutoff<1e-50 | FC cutoff > 2",
                gridlines.minor = FALSE,
                gridlines.major = FALSE,
                legendPosition = "bottom")
dev.off()



############ GSEA - Gene Set Enrichment Analysis
logger::log_info("Reading the Tabula Sapiens 1.0 skin file from the CellXgene repository")
#Select the genes with padj<0.01
dados%>%
  dplyr::filter(pvalue<0.01)->dados2
dados2$log2FoldChange->original_gene_list
dados2$entrez->names(original_gene_list)

# omit any NA values 
gene_list<-original_gene_list
gene_list<-gene_list[which(!is.na(names(gene_list)))]
gene_list<-gene_list[which(!duplicated(names(gene_list)))]


# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE) #5677 elements 5407 elements
set.seed(42) #Set seed for reproducibility
#Run GSEA
gse <- gseGO(geneList=gene_list, 
             ont ="BP", #"ALL"
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 270,  #Start with 1% of elements and set later
             maxGSSize = 800, #Start with 1000% of elements and set later
             pvalueCutoff = 0.1, 
             verbose = TRUE, 
             seed=TRUE,
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "BH") #"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"

logger::log_info("Plotting the GSEA results")
#Plotting the GSEA results 
png(file.path(fig_dir, "Fig1G_TBS_skin_Endo_by_NFATC1_exp_gse_5reg_endothelial subset_BP.png"),width=2800,height = 2500,res=600)
dotplot(gse, showCategory=5, split=".sign" ,color="pvalue") +
 facet_grid(.~.sign)+
 theme(text = element_text(size=15),
      panel.grid= element_blank(),
      axis.text.x = element_text(size=12,angle=90, hjust=1),
      axis.text.y = element_text(size=9,hjust=1))+xlim(c(0.35,0.55))
dev.off()
pdf(file.path(fig_dir, "Fig1G_TBS_skin_Endo_by_NFATC1_exp_gse_5reg_endothelial subset_BP.pdf"),width=4,height =5)
dotplot(gse, showCategory=5, split=".sign" ,color="pvalue") +
 facet_grid(.~.sign)+
 theme(text = element_text(size=15),
      panel.grid= element_blank(),
      axis.text.x = element_text(size=12,angle=90, hjust=1),
      axis.text.y = element_text(size=9,hjust=1))+xlim(c(0.35,0.55))
dev.off()


png(file.path(results_dir, "gse_20reg_endothelial subset_BP.png"),width=2500,height = 5200,res=600)
dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=12,angle=90, hjust=1),
        axis.text.y = element_text(size=9,hjust=1))

dev.off()

y <- setReadable(gse, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.csv(as.data.frame(y),file.path(results_dir,"gse.csv"), row.names=FALSE, quote=FALSE) 
#Plotting vascular development pathway rank
png(file.path(results_dir, "GSEA_sig6_select_pathway6_vdev.png"),width=4000,height = 3500,res=600)
gseaplot(gse, geneSetID = 4, title = gsub("^(.{27})(.*)$","\\1\n\\2",gse$Description[4]))
dev.off()


