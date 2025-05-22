
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

##Set work directory
setwd("your-path-here")


##Read file

##Data downloaded from Tabula Sapiens - Skin CellxGene repo
#(https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5)

readRDS("9405995f-b1c6-485c-8cc9-b80b9e870057.rds")->data_SKin

##Save first plot, dimplot
png("Resultados/DimPlot.png",width=7500,height=3000,res=600)
DimPlot(data_SKin, reduction = "umap",group.by = "cell_ontology_class")+
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

png("Resultados/Exp_cell.png",width=4500,height=3500,res=600)
DotPlot(object = data_SKin, features = "ENSG00000131196")+
  scale_x_discrete(labels=c( "NFATC1"))+xlab("")
dev.off()

#### Selecting endothelial cells

rownames(data_SKin@meta.data)[which(data_SKin@meta.data$cell_ontology_class=="endothelial cell")]

counts <- GetAssayData(data_SKin, slot="counts", assay="RNA") 

counts.sub <- counts[,which(data_SKin@meta.data$cell_ontology_class=="endothelial cell")] #only endothelial cells

metas.sub<-data_SKin@meta.data[which(data_SKin@meta.data$cell_ontology_class=="endothelial cell"),]

data_SKin2 <- CreateSeuratObject(counts=counts.sub,meta.data = metas.sub)

###Dichotomizing dataset
median(data_SKin2@assays$RNA$counts["ENSG00000131196",])->med_NFATC1

DimPlot(data_SKin2, reduction = "umap",group.by = "cell_ontology_class")+
  labs(title = "Tabula Sapiens (Skin)")

### Boxplot of log10 NFATC1 levels
boxplot(log10(data_SKin2@assays$RNA$counts["ENSG00000131196",]+1))


##Estimating cells expressing (or not) NFATC1
data.frame(NFATC1=data_SKin2@assays$RNA$counts["ENSG00000131196",])->exp

exp<-exp%>%
  mutate(Group=ifelse(NFATC1==0,"No expressing NFATC1","expressing NFATC1"))

table(exp$Group)

# expressing NFATC1 No expressing NFATC1 
# 510                  799 


exp$Group->data_SKin2@meta.data$groups


###Calculating DESEQ manually
rm(list=c("a","dados","exp","counts.sub","metas.sub","counts"))



####OTHER features
###Converting to a DESeq2 object

counts_ls <- list()
counts_ls <- data_SKin2@assays$RNA$counts+1
#names(counts_ls)<- "Tabula Skin"


str(counts_ls)

boxplot(log10(rowSums(counts_ls)))

log10(sum(counts_ls["ENSG00000131196",])) #4.403 NFATC1 sum level 

## In this case, no exclusion of genes will be applied
counts_ls[which(log10(rowSums(counts_ls))>3.15),]->counts_ls ## 3.11 is log10(1309) the number of participant cells


dds <- DESeqDataSetFromMatrix(counts_ls, 
                              colData=data_SKin2@meta.data,
                              design = ~groups)

### running differential expression analysis
dds <- DESeq(dds)
res<-results(dds, tidy=TRUE)
write.csv(res,"results_Tidy.csv") ##Save the differential expression result


############ Recall data


read.csv("results_Tidy.csv")->dados

##Converting Ensembl codes to gene symbols
library(EnsDb.Hsapiens.v79)

ensembl.genes <- dados$row

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID","ENTREZID"))

colnames(geneIDs1)<-c("Symbol","row","entrez")

merge(dados,geneIDs1,by="row")->dados

write.csv(dados[,-2],"Gene_List.csv") ##Saving diff expressed genes

############ Volcano plot

dados$log2FoldChange*-1 -> dados$log2FoldChange
png("Resultados/VolcanoPlot_Tabula_endothelial.png",width=4500,height = 3500,res=600)
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
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE) #5677 elements

gse <- gseGO(geneList=gene_list, 
             ont ="BP", #"ALL"
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 70,  #Start with 1% of elements and set later
             maxGSSize = 700, #Start with 1000% of elements and set later
             pvalueCutoff = 0.1, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "BH") #"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"

#View(as.data.frame(gse)) ##See results

png("Resultados/gse_5reg_endothelial subset_BP.png",width=2800,height = 2500,res=600)

dotplot(gse, showCategory=5, split=".sign" ,color="pvalue"
        ) + facet_grid(.~.sign)+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=12,angle=90, hjust=1),
        axis.text.y = element_text(size=9,hjust=1))+xlim(c(0.35,0.55))

dev.off()


png("Resultados/gse_20reg_endothelial subset_BP.png",width=2500,height = 5200,res=600)

dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=12,angle=90, hjust=1),
        axis.text.y = element_text(size=9,hjust=1))

dev.off()




y <- setReadable(gse, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.csv(as.data.frame(y),"gse.csv")

png("Resultados/select_pathway6.png",width=4000,height = 3500,res=600)
gseaplot(gse, geneSetID = 94, title = gsub("^(.{27})(.*)$","\\1\n\\2",gse$Description[94]))
dev.off()



