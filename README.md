# mcbarycenter

This package implements the marginal corrected barycenter (MCB) estimator.

At a high level, it estimates mixture distributions across thresholded outcomes
and combines them to estimate a quantile curve with uncertainty intervals.

Available mixture methods:

- `efron`
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
