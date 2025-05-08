# Analysis of **_NFATC1_** expression in Endothelial cells from single cell data - Tabula Sapiens

To analyze the expression levels of **_NFATC1_** in endothelial cells in the skin, we use the Tabula Sapiens dataset[](). This dataset contains single-cell RNA sequencing data from various tissues and cell types.


## Dependencies 

### Software

The following software is required to be installed and visible in the path before running the scripts:
- **R**: R `4.3.3`[**here**](https://cran.r-project.org/)

#### R packages and environment
To restore the R environment this can be done using the `renv` package. 

``` R
install.packages("renv")
library(renv)
renv::restore()

```
### Datasets used
All the analysis were perfomed using the Tabula Sapiens dataset. The following datasets were used:
- **Tabula Sapiens - endothelial cells from different tissues dataset**: This dataset contains single-cell RNA sequencing data from various tissues and cell types. It can be downloaded from the [Tabula Sapiens website](https://tabula-sapiens.sf.czbiohub.org/) with dataset ID `4676160f-b0f6-4dea-a104-c2ef86c74674`. The dataset can be accessed directly using the following URL: [4676160f-b0f6-4dea-a104-c2ef86c74674.rds](https://datasets.cellxgene.cziscience.com/4676160f-b0f6-4dea-a104-c2ef86c74674.rds).
- **Tabula Sapiens Skin dataset**: This dataset contains single-cell RNA sequencing data from skin tissue. It can be downloaded from the [Tabula Sapiens website](https://tabula-sapiens.sf.czbiohub.org/) with dataset ID `9405995f-b1c6-485c-8cc9-b80b9e870057`. The dataset can be accessed directly using the following URL: [9405995f-b1c6-485c-8cc9-b80b9e870057.rds](https://datasets.cellxgene.cziscience.com/9405995f-b1c6-485c-8cc9-b80b9e870057.rds).




