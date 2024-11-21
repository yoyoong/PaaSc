# PaaSc: Pathway Activity Analysis for Single Cells

## Introduction
***
The pathway is a network of interactions of multiple molecules in the cell, which is involved in regulating the 
physiological processes of the cell. By calculating the single-cell pathway activity score, we can evaluate the 
pathway activity level of a single cell. Here, we propose a new method for calculating pathway activity, PaaSc, 
based on multiple correspondence analysis and regression analysis. PaaSc can more accurately quantify the activity
of a pathway in a single cell and obtain the metabolic activity characteristics of the single-cell.

## Installation
PaaSc can be install with >R4.0 and support analysis under Seurat4 or Seurat5 (recommend).
***
Install PaaSc from github:
```
install.packages("devtools")
devtools::install_github("yoyoong/PaaSc", ref = "main")
```
Or install from local file:
```
download.file("https://github.com/yoyoong/PaaSc/releases/download/1.0.0/PaaSc_1.0.0.tar.gz","PaaSc_1.0.0.tar.gz") 
install.packages(PaaSc_1.0.0.tar.gz, repos = NULL, type="source")
```

## Usage
***
Two examples demonstrate how to use PaaSc to analyze scRNA-seq data and spatial RNA-seq data.

#### [Example1: pbmc3k scRNA-seq analysis](https://htmlpreview.github.io/?https://github.com/yoyoong/PaaSc/blob/main/example/Example1.html)
This example analyzes a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics.
#### [Example2: mouse brain spatial RNA-seq analysis](https://htmlpreview.github.io/?https://github.com/yoyoong/PaaSc/blob/main/example/Example2.html)
This example analyzes a spatial RNA-seq data of sagital mouse brain slices.

## References
***


## License
***
This project is licensed under the [GNU General Public License 3](LICENSE).