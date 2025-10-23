#' K-means Consensus Clustering
#'
#' Performs k-means consensus clustering using bootstrap resampling to generate
#' stable cluster assignments. This function runs multiple k-means clustering
#' iterations on bootstrapped samples of the data, builds a consensus matrix
#' based on co-clustering frequencies, and performs hierarchical clustering
#' on the consensus matrix.
#'
#' @param data A numeric matrix where rows are features and columns are samples.
#'   Row names are required for feature identification.
#' @param k Integer. The number of clusters to identify.
#' @param reps Integer. The number of bootstrap resampling iterations to perform.
#' @param ncores Integer. The number of CPU cores to use for parallel processing.
#' @param seed Integer. Random seed for reproducibility. Default is 0.
#'
#' @return An integer vector of cluster assignments for each feature, named with
#'   the feature names from the input data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Generates all unique feature pairs (including self-pairs)
#'   \item For each bootstrap iteration:
#'     \itemize{
#'       \item Samples columns (samples) with replacement
#'       \item Runs k-means clustering
#'       \item Records which feature pairs cluster together
#'     }
#'   \item Aggregates co-clustering frequencies across all iterations
#'   \item Builds a symmetric consensus matrix
#'   \item Performs hierarchical clustering (complete linkage) on 1 - consensus
#'   \item Cuts the dendrogram to obtain k final clusters
#' }
#'
#' The function uses parallel processing for both the bootstrap iterations and
#' the aggregation of results. Random number generation uses L'Ecuyer-CMRG
#' for reproducibility in parallel contexts.
#'
#' @importFrom dplyr slice_sample left_join select mutate if_else %>%
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom Matrix sparseVector
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ clusterSplit
#' @importFrom parallel nextRNGStream
#' @importFrom pbapply pblapply
#' @importFrom purrr reduce
#' @importFrom fastcluster hclust
#' @importFrom stats kmeans cutree as.dist
#' @importFrom methods as
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data matrix
#' set.seed(123)
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(data) <- paste0("Feature_", 1:100)
#'
#' # Run consensus clustering
#' clusters <- kmeansCC(
#'   data = data,
#'   k = 4,
#'   reps = 100,
#'   ncores = 2,
#'   seed = 0
#' )
#'
#' # View cluster assignments
#' table(clusters)
#' }
kmeansCC <- function(data, k, reps, ncores, seed = 0) {

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
    select(combination.1 = .data$key.x, combination.2 = .data$key.y)

  # Inner function: single bootstrap clustering iteration
  run_kmeans_iteration <- function(seed) {
    assign(".Random.seed", seed, envir = .GlobalEnv)

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
      dplyr::select(-.data$clustMember)

    # Determine co-clustering for all pairs
    out <- combinations2 %>%
      dplyr::left_join(clust_assignments, by = c("combination.1" = "key")) %>%
      dplyr::select(-.data$combination.1) %>%
      dplyr::left_join(clust_assignments, by = c("combination.2" = "key")) %>%
      dplyr::select(-.data$combination.2) %>%
      dplyr::mutate(coclustered = as.integer(dplyr::if_else(.[[1]] == .[[2]], 1, 0))) %>%
      dplyr::select(.data$coclustered)

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
    pivot_wider(names_from = .data$combination.2, values_from = .data$value) %>%
    column_to_rownames("combination.1")

  upper_coclust[lower.tri(upper_coclust, diag = FALSE)] <- 0
  lower_coclust <- t(upper_coclust)
  lower_coclust[upper.tri(lower_coclust, diag = TRUE)] <- 0
  coclust <- upper_coclust + lower_coclust

  # Hierarchical clustering on consensus matrix
  hc <- fastcluster::hclust(as.dist(1 - coclust), method = "complete")
  cutree(tree = hc, k = k)
}
