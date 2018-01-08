library("testthat")
library("knnsmoother")

context('knnsmoother: is output correct?')
num_genes <- 5
num_samples <- 10
x <- matrix(abs(sin(seq(from=1, to=num_samples*num_genes))) * 100,
            ncol = num_samples, byrow = F)

colnames(x) <- paste0('cell', seq_len(ncol(x)))
rownames(x) <- paste0('gene', seq_len(nrow(x)))

xd <- as.matrix(dist(t(x)))

source(file.path('src'))
