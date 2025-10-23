# Load the package function
library(ConsensusCluster)

(load("testData/RD8a_m1_imputed.RData"))


# z-score
exprs(m) <- t(scale(t(exprs(m))))
m1 <- m[1:100, 1:20]

# kmeans Consensus Clustering ---------------------------------------------

res <- kmeansCC(data = exprs(m1),
         k = 4,
         reps = 100,
         ncores = 4,
         seed = 0)
res
