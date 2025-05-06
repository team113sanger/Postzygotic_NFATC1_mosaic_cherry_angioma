#!/usr/bin/env Rscript
#Get the current working directory
here::i_am("scripts/tsapiens_analysis/tabula_sapiens_skin_analysis.R")
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
#Check if the results directory exists, if not create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
}

## Get the data from the Tabula Sapiens 
# Tabula Sapiens 1.0 Skin  URL https://datasets.cellxgene.cziscience.com/9405995f-b1c6-485c-8cc9-b80b9e870057.h5ad
# Download  .RDS with the Skin data from Tabula Sapiens 1.0 - ID 9405995f-b1c6-485c-8cc9-b80b9e870057.rds
# hd5ad files generated for TBS1.0 have issues in being converted to h5Seurat format to be loaded in Seurat
# Because of this, we are using the .rds file directly from the CellXgene repository
tbs_skin_data_url<-"https://datasets.cellxgene.cziscience.com/9405995f-b1c6-485c-8cc9-b80b9e870057.rds"
logger::log_info("Reading the Tabula Sapiens 1.0 skin file from the CellXgene repository")
logger::log_info(paste0("Downloading : ", tbs_skin_data_url))
data_SKin <- readRDS(url(tbs_skin_data_url))
##Data downloaded from Tabula Sapiens - Skin CellxGene repo
#(https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5)
#tab_skin_data<- file.path(tabdatadir, "9405995f-b1c6-485c-8cc9-b80b9e870057.rds")
#data_SKin <- readRDS(tab_skin_data)

##Save first plot, dimplot
png(file.path(results_dir, "DimPlot.png"),width=7500,height=3000,res=600)
DimPlot(data_SKin, reduction = "umap",group.by = "cell_ontology_class")+
  labs(title = "Tabula Sapiens (Skin)")
dev.off()

# With Colours more friendly to colorblind people
#FBE183 #F7D144 #F4C00D #FAAA05 #F78C09 #E35D2A #C8403D mast cell
 #A7373F #B14054 #D55474 #E26781 #E77988 #E78E97 #E3A0A5
 #C38AA3 #A8759F #A26296 #924E8A #723B79 #4F447E #2A6494
 #237E91 #299783 #5AAB6B #92C051
ident_colors <- colorRampPalette(paletteer::paletteer_d("MetBrewer::Signac",14 ))(length(unique(data_SKin$cell_ontology_class))) #Get colors for each cell type
ident_colors[11]<-"#C8403D" # Reassign Endothelial cell color
ident_colors[7]<-"#E26781" # Reassign Mast cell color
names(ident_colors) <- unique(data_SKin$cell_ontology_class) #Set names for each color
png(file.path(results_dir, "DimPlot_colblind.png"),width=7500,height=3000,res=600)
DimPlot(data_SKin, reduction = "umap",group.by = "cell_ontology_class", cols=ident_colors)+
  labs(title = "Tabula Sapiens (Skin)")
dev.off()

# Grey everything Else Endothelial cell only
gidentcol<-rep("lightgrey",length(unique(data_SKin$cell_ontology_class)))
gidentcol[11]<-"#C8403D" # Re-assign Endothelial cell color
names(gidentcol) <- unique(data_SKin$cell_ontology_class) #Set names for each color
png(file.path(results_dir, "DimPlot_gendothelial.png"),width=7500,height=3000,res=600)
DimPlot(data_SKin, reduction = "umap",group.by = "cell_ontology_class", cols=gidentcol)+
  labs(title = "Tabula Sapiens (Skin)")
dev.off()

##Set identity
data_SKin <- SetIdent(data_SKin, value = "cell_ontology_class")

##Sorting data by NFATC1
a<-DotPlot(object = data_SKin, features = "ENSG00000131196")
a$data%>%arrange(-avg.exp)%>%select(id)%>%as.matrix()

data_SKin@active.ident <- factor(data_SKin@active.ident, 
                            levels=a$data%>%arrange(pct.exp)%>%select(id)%>%as.matrix()
)

png(file.path(results_dir, "NFATC1_Exp_TS_skin_cells.png"),width=4500,height=3500,res=600)
DotPlot(object = data_SKin, features = "ENSG00000131196")+
  scale_x_discrete(labels=c( "NFATC1"))+xlab("")
dev.off()
#plot with NFATC1 and VEGFR2 VEGFR1
png(file.path(results_dir, "NFATC1_VGFR1_2_Exp_TS_skin_cells.png"),width=4500,height=3500,res=600)
DotPlot(object = data_SKin, features = c("ENSG00000131196", "ENSG00000102755", "ENSG00000128052" ))+
  scale_x_discrete(labels=c( "NFATC1", "VEGFR1", "VEGFR2"))+xlab("")+theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

#### Selecting endothelial cells
rownames(data_SKin@meta.data)[which(data_SKin@meta.data$cell_ontology_class=="endothelial cell")]

counts <- GetAssayData(data_SKin, slot="counts", assay="RNA") 

counts.sub <- counts[,which(data_SKin@meta.data$cell_ontology_class=="endothelial cell")] #only endothelial cells

metas.sub<-data_SKin@meta.data[which(data_SKin@meta.data$cell_ontology_class=="endothelial cell"),]

data_SKin2 <- CreateSeuratObject(counts=counts.sub,meta.data = metas.sub)

###Dichotomizing dataset
med_NFATC1<- median(data_SKin2@assays$RNA$counts["ENSG00000131196",])

### Boxplot of log10 NFATC1 levels
png(file.path(results_dir, "Boxplot_NFATC1.png"),width=4500,height=3500,res=600)
boxplot(log10(data_SKin2@assays$RNA$counts["ENSG00000131196",]+1))
dev.off()

##Estimating cells expressing (or not) NFATC1
exp<-data.frame(NFATC1=data_SKin2@assays$RNA$counts["ENSG00000131196",])
exp<-exp%>%
  mutate(Group=ifelse(NFATC1==0,"Endothelial cell NFATC1neg","Endothelial cell NFATC1+"))
table(exp$Group)
# expressing NFATC1 No expressing NFATC1 
# 510                  799 
exp$Group->data_SKin2@meta.data$groups

#Subset endothelial cells only
data_Skin3<-subset(data_SKin, cells=rownames(data_SKin@meta.data)[which(data_SKin@meta.data$cell_ontology_class=="endothelial cell")])
data_Skin3@meta.data$groups<-exp$Group
data_Skin3<- SetIdent(data_Skin3, value = "groups")
#Plot  the expression NFATC1, VEGFR1, VEGFR2FeaturePLot with NFATC1 and VEGFR2 VEGFR1
png(file.path(results_dir, "Fplot_Endocells_NFATC1_VEGFR.png"),width=4500,height=3500,res=600)
FeaturePlot(data_Skin3, features = c("ENSG00000131196", "ENSG00000102755", "ENSG00000128052" ))
dev.off()

#Plot  the expression NFATC1, VEGFR1, VEGFR2 and MAPK1 (ERK) FeaturePLot with NFATC1 and VEGFR2 VEGFR1
# Relevant for proliferation
# https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2020.599281/full
png(file.path(results_dir, "Exp_cell_endo_NFATC_explblgroup.png"),width=4500,height=3500,res=600)
DotPlot(object = data_Skin3, features = c("ENSG00000131196", "ENSG00000102755", "ENSG00000128052", "ENSG00000100030" ))+
  scale_x_discrete(labels=c( "NFATC1", "VEGFR1", "VEGFR2", "MAPK1"))+xlab("")+theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()
png(file.path(results_dir, "Endocells_DimPlot_by_NFATC1_expgroup.png"),width=6500,height=3000,res=600)
DimPlot(data_Skin3, reduction = "umap", group.by = "groups")+
  labs(title = "Tabula Sapiens (Skin)")
dev.off()

###Calculating DESEQ manually
rm(list=c("a","dados","exp","counts.sub","metas.sub","counts"))

####OTHER features
###Converting to a DESeq2 object
counts_ls <- list()
counts_ls <- data_SKin2@assays$RNA$counts+1
#names(counts_ls)<- "Tabula Skin"

str(counts_ls)

png(file.path(results_dir, "Boxplot_counts_ls.png"),width=4500,height=3500,res=600)
boxplot(log10(rowSums(counts_ls)))
dev.off()

log10(sum(counts_ls["ENSG00000131196",])) #4.403 NFATC1 sum level 

## In this case, no exclusion of genes will be applied
counts_ls[which(log10(rowSums(counts_ls))>3.15),]->counts_ls ## 3.11 is log10(1309) the number of participant cells


dds <- DESeqDataSetFromMatrix(counts_ls, 
                              colData=data_SKin2@meta.data,
                              design = ~groups)

### running differential expression analysis
dds <- DESeq(dds)
res<-results(dds, tidy=TRUE)
write.csv(res,file.path(results_dir,"results_Tidy.csv")) ##Save the differential expression result


############ Recall data
dados<-read.csv(file.path(results_dir,"results_Tidy.csv"))

##Converting Ensembl codes to gene symbols
library(EnsDb.Hsapiens.v79)

ensembl.genes <- dados$row

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID","ENTREZID"))

colnames(geneIDs1)<-c("Symbol","row","entrez")

merge(dados,geneIDs1,by="row")->dados

write.csv(dados[,-2],file.path(results_dir, "Gene_List.csv")) ##Saving diff expressed genes

############ Volcano plot
dados$log2FoldChange*-1 -> dados$log2FoldChange
png(file.path(results_dir, "VolcanoPlot_Tabula_endothelial.png"),width=4500,height = 3500,res=600)
EnhancedVolcano(dados,lab=dados$Symbol,y="padj",x="log2FoldChange",
                pCutoff = 1e-50,
                drawConnectors = T,
                title = "NFATC1+ vs. NFATC1- (subset: endothelial cells)",
                subtitle = "p cutoff<1e-50 | FC cutoff > 2",
                legendPosition = "bottom")
dev.off()


############ GSEA

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

#View(as.data.frame(gse)) ##See results
png(file.path(results_dir, "gse_5reg_endothelial subset_BP.png"),width=2800,height = 2500,res=600)
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

write.csv(as.data.frame(y),file.path(results_dir,"gse.csv"))

png(file.path(results_dir, "select_pathway6.png"),width=4000,height = 3500,res=600)
gseaplot(gse, geneSetID = 94, title = gsub("^(.{27})(.*)$","\\1\n\\2",gse$Description[94]))
dev.off()

### Generate the plot of NFATC1 expression on endothelial cells across all tissues on the Tabula Sapiens v1.0
all_tissues_tbs_url<-"https://datasets.cellxgene.cziscience.com/4676160f-b0f6-4dea-a104-c2ef86c74674.rds"
logger::log_info("Reading the Tabula Sapiens 1.0 all organs RDS from the CellXgene repository")
tbs_v1_organs<-readRDS(url(all_tissues_tbs_url))

##Set identity of cells based on the tissue they belong to on the publication
tbs_v1_organs <- SetIdent(tbs_v1_organs, value = "tissue_in_publication")

##Set identity
a<-DotPlot(object = tbs_v1_organs, features = "ENSG00000131196")
a$data%>%arrange(-avg.exp)%>%select(id)%>%as.matrix()
tbs_v1_organs@active.ident <- factor(tbs_v1_organs@active.ident, 
                            levels=a$data%>%arrange(pct.exp)%>%select(id)%>%as.matrix()
)

##Sorting data by NFATC1
png(file.path(results_dir, "Exp_NFATC1_tissues_tbs1.png"),width=3000,height = 3500,res=600)
DotPlot(object = tbs_v1_organs, features = "ENSG00000131196")+
  scale_x_discrete(labels=c( "NFATC1"))+xlab("")
dev.off()
rm(tbs_v1_organs)
