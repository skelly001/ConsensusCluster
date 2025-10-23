library(tidyverse)
library(MSnSet.utils)
library(pbapply)
library(parallel)

load("testData/shotgun_topdown_int_modann_cnt_wres.RData")

fData(m)$feature <- rownames(fData(m))

# Check data quality
plotNA(m)
boxplot(exprs(m), outline = FALSE)
hist(rowMeans(exprs(m), na.rm = TRUE), 100)

# Impute missing values
if (!file.exists("testData/RD8a_m1_imputed.RData")) {
  exprs(m) <- t(completeObs(pca(as(m, "ExpressionSet"),
                                method = "svdImpute",
                                nPcs = min(dim(m)),
                                center = TRUE)))
  save(m, file = "testData/RD8a_m1_imputed.RData")
} else {
  load("testData/RD8a_m1_imputed.RData")
}

# Z-score normalization
exprs(m) <- t(scale(t(exprs(m))))
m1 <- m[1:100, 1:20]


# K-means Consensus Clustering --------------------------------------------
data = exprs(m1)
  k = 4
  reps = 100
  ncores = 4
  seed = 0

kmeansCC <- function(data, k, reps, ncores, seed) {

  # Generate all unique pairs including self-pairs
  n_features <- length(rownames(data))
  feature_names <- rownames(data)
  idx <- which(upper.tri(matrix(0, n_features, n_features), diag = TRUE),
               arr.ind = TRUE)

  combinations <- data.frame(
    combination.1 = feature_names[idx[, 1]],
    combination.2 = feature_names[idx[, 2]],
    stringsAsFactors = FALSE
  )

  # Map feature names to integer keys for faster joins
  key <- data.frame(
    feature = rownames(data),
    key = seq_len(nrow(data))
  )

  combinations2 <- combinations %>%
    left_join(key, by = c("combination.1" = "feature")) %>%
    left_join(key, by = c("combination.2" = "feature")) %>%
    select(combination.1 = key.x, combination.2 = key.y)

  # Inner function: single bootstrap clustering iteration
  run_kmeans_iteration <- function(seed) {
    .Random.seed <- seed

    # Bootstrap sample
    data_boot <- dplyr::slice_sample(as.data.frame(t(data)),
                              prop = 1,
                              replace = TRUE) %>%
      as.matrix() %>%
      t()

    # Cluster assignment
    clust_assignments <- stats::kmeans(data_boot,
                                       centers = k,
                                       nstart = 1,
                                       iter.max = 100)$cluster

    clust_assignments <- data.frame(
      clustMember = names(clust_assignments),
      clustID = as.integer(clust_assignments),
      row.names = NULL
    ) %>%
      dplyr::left_join(key, by = c("clustMember" = "feature")) %>%
      dplyr::select(-clustMember)

    # Determine co-clustering for all pairs
    out <- combinations2 %>%
      dplyr::left_join(clust_assignments, by = c("combination.1" = "key")) %>%
      dplyr::select(-combination.1) %>%
      dplyr::left_join(clust_assignments, by = c("combination.2" = "key")) %>%
      dplyr::select(-combination.2) %>%
      dplyr::mutate(coclustered = as.integer(dplyr::if_else(.[[1]] == .[[2]], 1, 0))) %>%
      dplyr::select(coclustered)

    as(out[[1]], "sparseVector")
  }

  # Generate reproducible seeds
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  seeds <- vector("list", reps)
  seeds[[1]] <- .Random.seed
  for (i in seq(2, reps)) {
    seeds[[i]] <- nextRNGStream(seeds[[i - 1]])
  }

  # Parallel bootstrap iterations
  multiproc_cl <- makeCluster(ncores)
  clusterExport(multiproc_cl,
                varlist = c("data", "k", "combinations2", "key", "%>%"),
                envir = environment())
  clusterEvalQ(multiproc_cl,
               library(Matrix, include.only = "sparseVector",
                       quietly = TRUE, verbose = FALSE))

  res_all <- pblapply(cl = multiproc_cl, X = seeds, FUN = run_kmeans_iteration)
  tryCatch(stopCluster(multiproc_cl))

  # Parallel summation of results
  multiproc_cl <- makeCluster(ncores)
  res_all <- clusterSplit(cl = multiproc_cl, seq = res_all)
  res_all <- pblapply(
    cl = multiproc_cl,
    X = res_all,
    FUN = function(x) {
      purrr::reduce(x, function(z, zz) {
        `+`(as(z, "integer"), as(zz, "integer"))
      })
    }
  )
  tryCatch(stopCluster(multiproc_cl))

  # Calculate consensus proportions
  res_all <- res_all %>%
    purrr::reduce(function(z, zz) {
      `+`(as(z, "integer"), as(zz, "integer"))
    }) %>%
    `/`(reps)

  # Build symmetric co-clustering matrix
  upper_coclust <- combinations %>%
    mutate(value = res_all) %>%
    pivot_wider(names_from = combination.2, values_from = value) %>%
    column_to_rownames("combination.1")

  upper_coclust[lower.tri(upper_coclust, diag = FALSE)] <- 0
  lower_coclust <- t(upper_coclust)
  lower_coclust[upper.tri(lower_coclust, diag = TRUE)] <- 0
  coclust <- upper_coclust + lower_coclust

  # Hierarchical clustering on consensus matrix
  hc <- fastcluster::hclust(as.dist(1 - coclust), method = "complete")
  cutree(tree = hc, k = k)
}


# Run consensus clustering
t <- Sys.time()

res <- kmeansCC(
  data = exprs(m1),
  k = 4,
  reps = 100,
  ncores = 4,
  seed = 0
)

Sys.time() - t
