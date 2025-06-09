# M4R-project-codes

## R Code Overview

The following R scripts are provided for various simulation studies and analyses described in the thesis:

### Standard Parameter Simulation Study (Chapter 4.1)
- **File**: `out_sample_y.R`
  - Description: R code for the standard parameters simulation study.

### Varying Parameters Simulation Studies (Chapters 4.2 - 4.5)
- **Files**:
  - `out_sample_features.R`: R code for varying features simulation study.
  - `out_sample_correlation.R`: R code for varying correlation simulation study.
  - `out_sample_sparsity.R`: R code for varying sparsity simulation study.
  - `out_sample_strength.R`: R code for varying signal strength simulation study.

### Comparison of Standard and Varying Parameter Studies
- **Files**:
  - `elastic_net.R`: R code for comparison with standard parameters.
  - `elastic_net_features.R`: R code for comparison with varying features.
  - `elastic_net_correlation.R`: R code for comparison with varying correlation.
  - `elastic_net_sparsity.R`: R code for comparison with varying sparsity.
  - `elastic_net_signal_strength.R`: R code for comparison with varying signal strength.

### Real Data Studies (Chapter 5)
- **Files**:
  - `real_data.R`: R code for real data analysis.
  - `elastic_net_real_data.R`: R code for elastic net analysis on real data.

## Clustering Methods
The following functions and packages are utilized for clustering methods in the analysis:

- **SVD Decomposition**: Built-in function `svd()`
- **QR Decomposition**: Built-in function `qr()`
- **K-Means Clustering**: `kmeans()` function from the `cluster` package
- **NMF (Non-negative Matrix Factorization)**: `nmf()` function from the `RcppML` package

## Model Fitting Functions
- **Lasso, Group Lasso, and Sparse Group Lasso**: `dfr_sgl.cv()` function from the `dfr` package
- **Elastic Net**: `cv.glmnet` function from the `glmnet` package
