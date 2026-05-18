# mcbarycenter

This package implements the marginal corrected barycenter (MCB) estimator.

## Install

```r
devtools::install_github("jpspeng/mcbarycenter")
```

## Using the package

Available mixture methods:

- `spline`
- `npmle`
- `raw`
- `beta`

Example:

```r
library(mcbarycenter)
data("sparse_dist_data")

mcbary(sparse_dist_data,
    id_col = "id",
    val_col = "value",
    method = "npmle",
    cutpoints = 20,
    bootstrap_samples = 20,
    use_midpoint = T,
    estimate_endpoints = F,
    ci_method = "wald")
```
