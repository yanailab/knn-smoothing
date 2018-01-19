# K-nearest neighbor smoothing for UMI-filtered scRNA-Seq data
# (R implementation)

# Author: Yun Yan <yun.yan@nyumc.org>
# Copyright (c) 2017 New York University

library(Matrix)

freeman_tukey_transform <- function(mat){
  sqrt(mat) + sqrt(mat + 1)
}

pdist <- function(mat){
  # @param mat A non-negative matrix with samples by features
  # @reference http://r.789695.n4.nabble.com/dist-function-in-R-is-very-slow-td4738317.html
  mtm <- Matrix::tcrossprod(mat)
  sq <- rowSums(mat^2)
  out0 <- outer(sq, sq, "+") - 2 * mtm
  out0[out0 < 0] <- 0

  sqrt(out0)
}

smoother_calculate_distances <- function(mat){
  # @param mat A matrix with genes by samples

  # normalize to median transcript count
  num_transcripts <- Matrix::colSums(mat)
  size_factor <- median(num_transcripts, na.rm = T) / num_transcripts

  mat_norm <- t(t(mat) * size_factor)
  # apply freeman-tukey transform
  mat_FTT <- freeman_tukey_transform(mat_norm)
  # calculate all pairwise distances using the Euclidean metric
  mat_D <- pdist(t(mat_FTT))

  as.matrix(mat_D)
}

smoother_aggregate_nearest_nb <- function(mat, D=NULL, k=5){
  # @param mat A matrix in a shape of #genes x #samples.
  # @param D A predefined distance matrix in a shape of #samples x #samples. If
  #   not specified, D is the distance matrix of the input \code{mat}.
  # @param k An integer to choose \code{k} nearest samples (self-inclusive) to
  #  aggregate based on the distance matrix \code{D}. If \code{k} is greater than
  #  #samples, \code{k} is forced to be #samples to continue aggregation.
  if (is.null(D)) D <- pdist(mat)

  sapply(seq_len(ncol(mat)), function(cid){
    nb_cid <- head(order(D[cid, ]), k)
    closest_mat <- mat[, nb_cid]
    return(Matrix::rowSums(closest_mat))
  })
}

#' KNN-smoothing on UMI-filtered single-cell RNA-seq data
#'
#' @param mat A numeric matrix with gene names on rows and cell names on columns.
#' @param k Number of nearest neighbours to aggregate.
#' @return A smoothed numeric matrix.
#' @examples
#' X <- matrix(abs(sin(seq(from=1, to=1000, length.out = 1000))),
#' nrow = 25, byrow = T)
#' y <- rep(1:4, each=10)
#' dim(X)
#' colnames(X) <- as.character(paste0("s", seq_len(ncol(X))))
#' rownames(X) <- as.character(paste0("g", seq_len(nrow(X))))
#' S <- knn_smoother(X, k=5)
#' plot(X[1, ], X[3, ], col=factor(y), main="original")
#' plot(S[1, ], S[3, ], col=factor(y), main="smoothed")
#' @export
knn_smoother <- function(mat, k=5){
  if (k >= ncol(mat)) stop('k should be less than the number of available samples.')

  cname <- colnames(mat)
  gname <- rownames(mat)

  num_powers <- ceiling(log2(k + 1))
  S <- mat
  for (p in seq(1, num_powers)){
    k_step <- min(2^p - 1, k)
    message(paste0('Step ', p, '/', num_powers, ':',
                   'Smoothing using k=', k_step))
    D <- smoother_calculate_distances(S)
    S <- smoother_aggregate_nearest_nb(mat, D, k_step + 1)
  }
  if (! is.null(cname)) colnames(S) <- cname
  if (! is.null(gname)) rownames(S) <- gname

  S
}
