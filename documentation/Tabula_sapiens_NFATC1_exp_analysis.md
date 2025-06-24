# Analysis of **_NFATC1_** expression in Endothelial cells from single cell data - Tabula Sapiens

To analyze the expression levels of **_NFATC1_** in endothelial cells within the skin and other tissues we use the Tabula Sapiens dataset. This dataset contains single-cell RNA sequencing data from various tissues and cell types.

## Dependencies 

### Software

The following software is required to be installed and visible in the path before running the scripts:
- **R**: R `4.3.3`[**here**](https://cran.r-project.org/)

#### R packages and environment
To restore the R environment used for the analysis this can be done using the `renv` package. If you're interested in [reproducing the R environment using `renv`](https://rstudio.github.io/renv/reference/index.html) please follow the official documentation. All the dependencies and versions used are listed in the `renv.lock` file. 

Example of how to restore the R environment using `renv`:
``` R
install.packages("renv")
library(renv)
renv::restore()
```

### Datasets used
All the analysis were performed using the Tabula Sapiens dataset. The following datasets were used:
- **Tabula Sapiens - endothelial cells from different tissues dataset**: This dataset contains single-cell RNA sequencing data from various tissues and cell types. It can be downloaded from the [Tabula Sapiens website](https://tabula-sapiens.sf.czbiohub.org/) with dataset ID `4676160f-b0f6-4dea-a104-c2ef86c74674`. The dataset can be accessed directly using the following URL: [4676160f-b0f6-4dea-a104-c2ef86c74674.rds](https://datasets.cellxgene.cziscience.com/4676160f-b0f6-4dea-a104-c2ef86c74674.rds).
- **Tabula Sapiens Skin dataset**: This dataset contains single-cell RNA sequencing data from skin tissue. It can be downloaded from the [Tabula Sapiens website](https://tabula-sapiens.sf.czbiohub.org/) with dataset ID `9405995f-b1c6-485c-8cc9-b80b9e870057`. The dataset can be accessed directly using the following URL: [9405995f-b1c6-485c-8cc9-b80b9e870057.rds](https://datasets.cellxgene.cziscience.com/9405995f-b1c6-485c-8cc9-b80b9e870057.rds).

## Analysis and figures reproduction

To run the analysis of **_NFATC1_** expression in endothelial cells run the following R script:

``` bash
Rscript 02_tabula_sapiens_endothelial_cells_analysis.R
```
This script will download the required datasets, perform _NFATC1_ expression analysis, and generate the figures. The results will be saved in the `results/figures` directory.

## Results

### Figures
The analysis will generate the figures in both PNG and PDF formats. The figures generated will include:

- **Fig1C_TBS_skin_DimPlot_colblind.pdf**: Containing the results of Fig1 C panel, which contains reduction plot of skin cells from Tabula Sapiens, coloured by cell type.
- **Fig1D_NFATC1_Exp_TS_skin_cells.pdf**: Containing the results of Fig1 D panel, which shows the expression of _NFATC1_ in skin cells from Tabula Sapiens.
- **Fig1E_Exp_NFATC1_tissues_tbs1.pdf**: Containing the results of Fig1 E panel, which shows the expression of _NFATC1_ in endothelial cells across different tissues from Tabula Sapiens.
- **Fig1F_VolcanoPlot_Tabula_endothelial.pdf**: Containing the results of Fig1 F panel, which shows the volcano plot of differentially expressed genes between _NFATC1_ positive and negative endothelial cells.
- **Fig1G_TBS_skin_Endo_by_NFATC1_exp_gse_5reg_endothelial subset_BP.pdf**: Containing the results of Fig1 G panel, which shows the set of Biological Processes (BP) enriched in the differentially expressed genes found between _NFATC1_ positive and endothelial cells from Tabula Sapiens skin dataset.
- **Fig1H_TBS_skin_Endothelial_cells_by_NFATC_expgroup.pdf**: Containing the results of Fig1 H panel, which shows the expression of _NFATC1_ in endothelial cells from Tabula Sapiens skin dataset, grouped by expression levels.

### Additional results
Additional plots and tables with the list of differentially expressed genes and their associated statistics will be saved in the `results/tsapiens_analysis` directory. 

```bash
tree results/tsapiens_analysis
results/tsapiens_analysis
├── Boxplot_counts_ls.png
├── Boxplot_NFATC1.png
├── DE_genes_by_NFATC1exp_results_Tidy_NFATC1_VEGFR1_2_MAPK1.tsv
├── DE_genes_by_NFATC1exp_results_Tidy.csv
├── gse_20reg_endothelial subset_BP.png
├── gse.csv
├── GSEA_sig6_select_pathway6_vdev.png
├── TB_skin_NFATC1_VGFR1_2_MAPK1_Exp_TS_skin_cells.png
├── TBS_skin_DimPlot_gendothelial.png
├── TBS_skin_Endocells_DimPlot_by_NFATC1_expgroup.png
└── TS_WB_endo_Exp_NFATC1_VEGFR1_2_MAPK1_tissues_tbs1.png
1 directory, 11 files
```