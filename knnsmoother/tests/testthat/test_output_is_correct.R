library("testthat")
library("knnsmoother")

context('knnsmoother: is output correct?')
e <- new.env() ## without messing up user's global env

e$num_genes <- 5
e$num_samples <- 10
e$x <- matrix(abs(sin(seq(from=1,
                          to=e$num_samples * e$num_genes))) * 100,
              ncol = e$num_samples, byrow = F)

colnames(e$x) <- paste0('cell', seq_len(ncol(e$x)))
rownames(e$x) <- paste0('gene', seq_len(nrow(e$x)))

e$xd <- as.matrix(dist(t(e$x)))

source(file.path("..", "..", "R", "essence.R"), local=e)

test_that("k=0: expect no smoothing effect at all.", {
  e$k <- 0
  expect_equal(c(e$x), c(knnsmoother::knn_smoothing(e$x, e$k)))
})

test_that("k=1: smoothing with only 1 neighbor.", {
  e$k <- 1
  expect_equal(c(e$r_knn_smoothing(e$x, e$k)),
               c(knnsmoother::knn_smoothing(e$x, e$k)))
})

test_that("k= half of #samples.",{
  e$k <- e$num_samples / 2
  expect_equal(c(e$r_knn_smoothing(e$x, e$k)),
               c(knnsmoother::knn_smoothing(e$x, e$k)))
})

test_that("k=#samples: illegal situation but it is considered full-set smoothing.",{
  e$k <- e$num_samples
  expect_equal(c(e$r_knn_smoothing(e$x, e$k)),
               c(knnsmoother::knn_smoothing(e$x, e$k)))
})

test_that("k= #samples + 1: illegal situation and cannot be allowed.",{
  e$k <- e$num_samples + 1
  expect_error(knnsmoother::knn_smoothing(e$x, e$k))
})
