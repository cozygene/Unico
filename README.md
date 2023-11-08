# Unico
There are the code implementation of a Unico, a unified cross-omics method designed to deconvolve standard 2-dimensional bulk matrices of samples by features into a 3-dimensional tensors representing samples by features by cell types. Unico stands out as the first principled model-based deconvolution method that is theoretically justified for any heterogeneous genomic data.

A R pacakge will be available soon on CRAN

# Installation Instructions
Please install the following packages in R, which should finish within an hour.
```
install.packages("pbapply")
install.packages("config")
install.packages("data.table")
install.packages("matrixStats")
install.packages("matrixcalc")
install.packages("mgcv")
install.packages("nloptr")
install.packages("testit")
install.packages("compositions")
install.packages("MASS")
install.packages("pracma")
```
# Version info
Our algorithm is tested on both R 3.6.1 and R 4.1.0, on both Linus and MacOS based machines.
The entirety of the system environment is also included at the end in the tutorial notebook and its associated html exported version.


# Running Unico
```
Please provide (1) Bulk matrix of raw un-normalized counts or methylation beta values as input: X (features by samples),
               (2) The (estimated) sample specific cell-type proportions matrix: W (samples by cell types).
                   This matrix should be normalized such that proportions on each row/individual sum up to 1
               (3) Optionally, the model can also handle cell-type level covariates: C1 (samples by number of cell-type level covariates)
               (4) Optionally, the model can also handle tissue-level covariates: C2 (samples by number of tissue-level covariates)
We encourage users to provide all matrices with explicit and meaningful row and column names when possible to avoid confusion regarding the dimension of the output object.
```
### First step is always to performing parameter estimation.
#### Unico will learn cell-type specific means, variances and covariances along with effect sizes on all covariates using GMM
```
unico.mdl = list()
unico.mdl$params.hat <- Unico(X, W, C1, C2, parallel = TRUE)
```
### Once parameters are estimated we can further estimate cell-type level profiles: cell types by features by samples, denoted as Z.hat
```
unico.mdl$Z.hat = tensor(X, W = W, C1, C2, unico.mdl$params.hat)
```
### Alternatively, we can perform parametric association study (under normality assumption) without explicitly derive the underlying cell-type levels:
```
Unico.mdl$params.hat = add_C1_C2_pvals_parametric(X = X, Unico.mdl = Unico.mdl$params.hat, slot_name = "parametric")
```
### Since the framework does not have any distribution assumption, we can also perform non-parametric association study (p-values are asymptotically derived):
```
Unico.mdl$params.hat = add_C1_C2_pvals_asymptotic(X = X, Unico.mdl = Unico.mdl$params.hat, slot_name = "asymptotic")
```
# Tutorial
Please head to "Tutorial/Tutorial.ipynb"
