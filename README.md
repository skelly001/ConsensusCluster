# ConsensusCluster

An R package for performing k-means consensus clustering with bootstrap resampling.

## Installation

You can install the package directly from this repository using `devtools`:

``` r
# Install devtools if you haven't already
install.packages("devtools")

# Install ConsensusCluster
devtools::install_github("skelly001/ConsensusCluster")
```

Or install locally:

``` r
devtools::install("path/to/ConsensusCluster")
```

## Quick Start

``` r
library(ConsensusCluster)

# Create sample data
set.seed(123)
data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(data_matrix) <- paste0("Feature", 1:100)

# Run consensus clustering
clusters <- kmeansCC(
  data = data_matrix,
  k = 4,
  reps = 100,
  ncores = 4,
  seed = 0
)

# View cluster assignments
table(clusters)
```

## Function: kmeansCC

Performs consensus clustering using k-means algorithm with bootstrap resampling.

### Parameters

-   `data`: A numeric matrix or data frame with features in rows and samples in columns
-   `k`: Integer. The number of clusters to identify
-   `reps`: Integer. The number of bootstrap iterations to perform (default: 100)
-   `ncores`: Integer. The number of CPU cores to use for parallel processing (default: 1)
-   `seed`: Integer. Random seed for reproducibility (default: 0)

### Returns

A named integer vector of cluster assignments for each feature.

### Algorithm

1.  Generate all unique pairs of features (including self-pairs)
2.  For each bootstrap iteration:
    -   Resample samples with replacement
    -   Perform k-means clustering
    -   Record which features cluster together
3.  Calculate consensus proportions: how often each pair of features co-clustered
4.  Build a symmetric co-clustering matrix
5.  Perform hierarchical clustering on 1 - consensus matrix
6.  Cut tree to obtain k final clusters

## License

MIT License