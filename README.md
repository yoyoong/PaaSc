# PaaSc: Pathway Activity Analysis for Single Cells

## Introduction
Cellular pathways represent complex networks of molecular interactions that regulate biological processes, from signal transduction to metabolic functions. The ability to accurately quantify pathway activity in individual cells within their spatial context represents a crucial step toward understanding cellular behavior, tissue organization, and disease mechanisms. Here, we propose a new method for calculating pathway activity, PaaSc, based on multiple correspondence analysis (MCA) and linear regression. PaaSc can more accurately quantify the activity of a pathway in scRNA-seq and spatial transcriptomics.

## Installation
PaaSc requires R version 4 or later and depends on Seurat (Version 4 and 5).

Install PaaSc from Github:
```
install.packages("devtools")
devtools::install_github("yoyoong/PaaSc", ref = "main")
```
Or install from local source:
```
download.file("https://github.com/yoyoong/PaaSc/releases/download/1.0.0/PaaSc_1.0.0.tar.gz","PaaSc_1.0.0.tar.gz") 
install.packages("PaaSc_1.0.0.tar.gz", repos = NULL, type="source")
```

## Tutorial

#### [Example 1: pbmc3k scRNA-seq](https://jiantaoshi.github.io/PaaSc/PaaSc_scRNAseq.html)
This tutorial demonstrates a typical workflow for analyzing pathway activities in scRNA-seq data using PaaSc.

#### [Example 2: mouse brain spatial RNA-seq](https://jiantaoshi.github.io/PaaSc/PaaSc_spatial.html)
This tutorial demonstrates a typical workflow for analyzing pathway activities in spatial RNA-seq data using PaaSc.

## References


## License
This project is licensed under the [GNU General Public License 3](LICENSE).
