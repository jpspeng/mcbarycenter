# mcbarycenter

This package implements the marginal corrected barycenter (MCB) estimator.

## Install

```{r}
devtools::install_github("jpspeng/npthreshold")
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
data("sample_data")

mcb(sample_data,
    id_col = "id",
    val_col = "value",
    method = "npmle",
    cutpoints = 20,
    bootstrap_samples = 20,
    use_midpoint = T,
    estimate_first_last = F,
    ci_method = "wald")
```
