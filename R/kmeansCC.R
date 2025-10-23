library(tidyverse)
library(MSnSet.utils)
library(glue)
library(magrittr)
# library(fastcluster)
# library(RcppAlgos)
library(pbapply)
library(parallel)


(load("testData/shotgun_topdown_int_modann_cnt_wres.RData"))

# add fData for making remaking rownames
fData(m)$feature <- rownames(fData(m))

# check data
plotNA(m) #  must impute
boxplot(exprs(m), outline = F)
rowMeans(exprs(m), na.rm = T) %>% hist(100)

# # Impute
if (!file.exists("testData/RD8a_m1_imputed.RData")) {
   exprs(m) <- t(pcaMethods::completeObs(pcaMethods::pca(as(m,"ExpressionSet"), method="svdImpute",
                                  nPcs=min(dim(m)), center = TRUE)))
   save(m, file = "testData/RD8a_m1_imputed.RData")
} else {
   (load("testData/RD8a_m1_imputed.RData"))
}


# z-score
exprs(m) <- t(scale(t(exprs(m))))


# kmeans Consensus Clustering ---------------------------------------------


data <- exprs(m1)
k <- 200
reps <- 10000
ncores <- 64
seed <-  0

t <- Sys.time()


combinations <- RcppAlgos::comboGeneral(rownames(data), 2, nThreads = 1, repetition = T) %>%
   as.data.frame()

colnames(combinations) <- c("combination.1", "combination.2")

######### must be fast
# combinations <- combn(rownames(data), 2) %>% 
#    t() %>% 
#    as.data.frame()
# 
# colnames(combinations) <- c("combination.1", "combination.2")
# 
# combinations <- combinations %>% 
#    bind_rows(data.frame("combination.1" = rownames(data),
#                         "combination.2" = rownames(data)))

########
# remove duplicate combinations (a, b == b, a) ; diag will be left
# combinations <- expand_grid(rownames(data), rownames(data))

# features <- data.frame(feat = rownames(data)) %>% 
#    arrange(feat)
# 
# ff <- function(x, features, len) { features[1:(len-x), ] }
# 
# o <- map(.x = 0:(length(rownames(data))-1),
#          .f = ff,
#          features = features,
#          len = length(rownames(data)))
# 
# # f <- \(x, idx, feat, len){
# #   filter(x, combination.2 %in% feat[1:(len-(idx-1))])
# # }
# f <- \(x, y, feat, len){
#    filter(x, combination.2 %in% y)
# }
# 
# combinations <- combinations %>% 
#    arrange("combination.1") %>% 
#    nest(data = "combination.2") %>% 
#    mutate(data = map2(.x = data,
#                       .y = o,
#                       .f = f,
#                       feat = features[[1]],
#                       len = length(rownames(data)))) %>% 
#    unnest(cols = data)
# colnames(combinations) <- c("combination.1", "combination.2")

###########
# combinations <- expand.grid(rownames(data), rownames(data))
# combinations <- combinations %>%
#    filter(Var1 != Var2)
# combinations <- combinations %>%
#    mutate(combo_id = case_when(
#       as.numeric(Var1) > as.numeric(Var2) ~ paste0(Var1, Var2),
#       .default = paste0(Var2, Var1)))
# combinations <- combinations[!duplicated(combinations$combo_id),]
# combinations <- combinations %>%
#    dplyr::select(-combo_id) %>%
#    bind_rows(data.frame(Var1 = rownames(data),
#                         Var2 = rownames(data)))
# colnames(combinations) <- c("combination.1", "combination.2")





key <- data.frame("feature" = rownames(data),
                  "key" = 1:length(rownames(data)))


combinations2 <- combinations %>% 
   left_join(key, by = c("combination.1" = "feature")) %>% 
   left_join(key, by = c("combination.2" = "feature")) %>% 
   select(colnames(.)[3:4])

colnames(combinations2) <- c("combination.1", "combination.2")



kmeansCC <- function(x) {
   
   # .Random.seed <- seed
   
   data_boot <- dplyr::slice_sample(as.data.frame(t(data)), prop = 1, 
                                    replace = T) %>% 
      as.matrix() %>% 
      t()
   
   clust_assignments <- stats::kmeans(data_boot,
                                      centers = k,
                                      nstart = 1,
                                      iter.max = 100)$cluster
   
   
   clust_assignments <- data.frame(clustMember = names(clust_assignments),
                                   clustID = as.integer(clust_assignments),
                                   row.names = NULL) %>%
      dplyr::left_join(key, by = c("clustMember" = "feature")) %>%
      dplyr::select(-clustMember)
   
   
   out <- combinations2 %>%
      dplyr::left_join(clust_assignments, by = c("combination.1" = "key")) %>%
      dplyr::select(-combination.1) %>%
      dplyr::left_join(clust_assignments, by = c("combination.2" = "key")) %>%
      dplyr::select(-combination.2) %>%
      dplyr::mutate(coclustered = as.integer(dplyr::if_else(.[1] == .[2], 1, 0))) %>%
      dplyr::select(coclustered)
   
   out <- as(out[[1]], "sparseVector")
   
   return(out)
   
}


multiproc_cl <- makeCluster(ncores)

# RNGkind("L'Ecuyer-CMRG")
# 
# set.seed(seed)
# 
# seeds <- vector("list", reps)
# seeds[[1]] <- .Random.seed
# 
# for (i in seq(from = 2, to = reps)) {
#    seeds[[i]] <- nextRNGStream(seeds[[i - 1]])
#    }

clusterSetRNGStream(multiproc_cl, iseed = seed)

clusterExport(multiproc_cl, varlist = c("data", "k", "combinations2", 
                                        "key", "%>%"),
              envir = environment())

clusterEvalQ(cl = multiproc_cl, expr = library(Matrix, 
                                               include.only = "sparseVector", 
                                               quietly = TRUE,
                                               verbose = FALSE))

res_all <- pblapply(cl = multiproc_cl, X = 1:reps, FUN = kmeansCC)


tryCatch(stopCluster(multiproc_cl))






multiproc_cl <- makeCluster(ncores)

res_all <- clusterSplit(cl = multiproc_cl, seq = res_all)

res_all <- pblapply(cl = multiproc_cl, X = res_all, 
                    FUN =  \(x){purrr::reduce(.x = x, 
                                              .f = \(z, zz){`+`(as(z, "integer"),
                                                                as(zz, "integer"))})})

tryCatch(stopCluster(multiproc_cl))

res_all <- res_all %>% 
   purrr::reduce(.f = \(z, zz){`+`(as(z, "integer"),
                                   as(zz, "integer"))}) %>% 
   `/`(reps)



upper_coclust <- combinations %>% 
   mutate(value = res_all) %>% 
   pivot_wider(names_from = combination.2,
               values_from = value) %>% 
   column_to_rownames("combination.1")
   


# lower tri NA to 0
upper_coclust[lower.tri(upper_coclust, diag = FALSE)] <- 0

# create lower matrix
lower_coclust <- t(upper_coclust)

# eliminate diag from lower
lower_coclust[upper.tri(lower_coclust, diag = TRUE)] <- 0


# symetric coclustering matrix
coclust <- upper_coclust + lower_coclust



hc <- fastcluster::hclust(d = as.dist( 1 - coclust ), method = "complete")

ct = cutree(tree = hc, k = k)
   
   

Sys.time() - t


