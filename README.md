# Unico
  <!-- badges: start -->
  [![R-CMD-check](https://github.com/cozygene/Unico/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cozygene/Unico/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

We present Unico, a unified cross-omics method designed to deconvolve standard 2-dimensional bulk matrices of samples by features into a 3-dimensional tensors representing samples by features by cell types. Unico stands out as the first principled model-based deconvolution method that is theoretically justified for any heterogeneous genomic data.

All necessary scripts used for the analyses reported in the manuscript can be found under folder Rscripts

Unico will be available soon on CRAN as an R package

### Install from GitHub
```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("https://github.com/cozygene/Unico")
```

### Version info
Our package is tested on both R 3.6.1 and R 4.1.0, on Windows, Linux and MacOS based machines.

### Tutorial
Please head to <a href="https://cozygene.github.io/Unico/articles/Unico-Tutorial.html">this vignette</a> for a step by step tutorial on (1) deconvolving a simulated PBMC pseudo-bulk expression dataset and (2) association testing on subsets of publicly available methylation datasets.

### Author

This software is developed by Zeyuan Johnson Chen (johnsonchen@cs.ucla.edu) and Elior Rahmani (EliorRahmani@mednet.ucla.edu).

### License

Unico is available under the <a href="https://opensource.org/license/gpl-3-0" target="_blank">GPL-3 license</a>.
