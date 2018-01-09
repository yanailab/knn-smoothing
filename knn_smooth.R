# K-nearest neighbor smoothing for UMI-filtered scRNA-Seq data
# (R implementation)

# Author: Yun Yan <yun.yan@nyumc.org>
# Copyright (c) 2017 New York University

library(Matrix)
library(magrittr)

r_freeman_tukey_transform <- function(mat){
  sqrt(mat) + sqrt(mat + 1)
}


r_calculate_distances <- function(mat){
  # mat: gene by sample
  # normalize to median transcript count
  num_transcripts <- Matrix::colSums(mat)
  size_factor <- median(num_transcripts, na.rm = T) / num_transcripts

  mat_norm <- t(t(mat) * size_factor)
  # apply freeman-tukey transform
  mat_FTT <- r_freeman_tukey_transform(mat_norm)
  # calculate all pairwise distances using the Euclidean metric
  mat_D <- dist(t(mat_FTT), method = "euclidean",
               upper = T, diag = T)
  return(as.matrix(mat_D))
}


#" KNN-smoothing on UMI-filtered single-cell RNA-seq data
#"
#" @param mat A numeric matrix with gene names on rows and cell names on columns.
#" @param k Number of nearest neighbours to aggregate.
#" @return A smoothed numeric matrix.
#" @examples
#" X <- matrix(abs(sin(seq(from=1, to=1000, length.out = 1000))),
#" nrow = 25, byrow = T)
#" y <- rep(1:4, each=10)
#" dim(X)
#" colnames(X) <- as.character(paste0("s", seq_len(ncol(X))))
#" rownames(X) <- as.character(paste0("g", seq_len(nrow(X))))
#" S <- knn_smoothing(X, k=5)
#" plot(X[1, ], X[3, ], col=factor(y), main="original")
#" plot(S[1, ], S[3, ], col=factor(y), main="smoothed")
r_knn_smoothing <- function(mat, k=5){
  if (k > ncol(mat)) stop('k should not be greater than the number of available samples')
  cname <- colnames(mat)
  gname <- rownames(mat)
  num_powers <- ceiling(log2(k + 1))
  S <- mat
  for (p in seq(1, num_powers)){
    k_step <- min(2^p - 1, k)
    message(paste0('Step ', p, '/', num_powers, ':',
                   'Smoothing using k=', k_step))
    D <- r_calculate_distances(S)
    S <- sapply(cname, function(cn){
      closest_id <- D[cn, ] %>% sort(.) %>%
        head(., k_step+1) %>%
        names(.)
      closest_mat <- mat[gname, closest_id] %>%
        matrix(., nrow=length(gname), byrow = F)
      rownames(closest_mat) <- gname
      colnames(closest_mat) <- closest_id
      return(Matrix::rowSums(closest_mat))
    })
  }
  S
}
