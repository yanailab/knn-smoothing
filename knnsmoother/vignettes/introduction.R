## ------------------------------------------------------------------------
library(knnsmoother)
num_genes <- 20
num_cells <- 100
x <- matrix(sample.int(num_genes*num_cells, replace = T), ncol = num_cells)
rownames(x) <- paste0('g', seq_len(nrow(x)))
colnames(x) <- paste0('c', seq_len(ncol(x)))
smoothed_x <- knn_smoothing(X = x, k = 5)

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(knnsmoother))
suppressPackageStartupMessages(library(microbenchmark))

## ------------------------------------------------------------------------
# pbmc4k
http_pbmc4k <- 'http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz'

## ------------------------------------------------------------------------
demo_dir <- tempdir()
savetotar <- file.path(demo_dir, 'demo_data.tar.gz')
cmd <- paste('wget -nv ', http_pbmc4k,
             '-O', savetotar)
system(cmd)

## ------------------------------------------------------------------------
demo_data_dir <- file.path(demo_dir, 'demo_data')
dir.create(demo_data_dir)
cmd <- paste('tar -zxvf', savetotar,
             '--directory', demo_data_dir)
system(cmd)

## ------------------------------------------------------------------------
# expr_spmat <- Read10X(file.path(
#   demo_data_dir, 'filtered_gene_bc_matrices', 'GRCh38'))
expr_spmat <- Matrix::readMM(file.path(
  demo_data_dir, 'filtered_gene_bc_matrices', 'GRCh38', 'matrix.mtx'))

## ------------------------------------------------------------------------
print(class(expr_spmat))
expr_mat <- as.matrix(expr_spmat); rm(expr_spmat);
print(dim(expr_mat))

## ------------------------------------------------------------------------
g_use_idx <- rowSums(expr_mat) > 0
c_use_idx <- colSums(expr_mat) >= 100
expr_mat <- expr_mat[g_use_idx, c_use_idx]
print(paste(nrow(expr_mat), 'genes', 'x', ncol(expr_mat), 'cells'))

