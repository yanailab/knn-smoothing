#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix cxxfreeman_tukey_transform(NumericMatrix X){
  // Input:
  // - X: a matrix in shape of gene by sample
  // Return:
  // a matrix in a same shape with freeman tukey transformation
  // (sqrt(X) + sqrt(X + 1))
  mat X_mat = as<mat>(X);
  return wrap(sqrt(X_mat) + sqrt(X_mat + 1));
}

// Equivalent to R's `dist(method = "euclidean")` function
// http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/
// [[Rcpp::export]]
NumericMatrix cxxdist_euclidean(NumericMatrix Ar, NumericMatrix Br) {
  // Input:
  // - Ar: a matrix with samples as rows
  // - Br: a matrix with samples as rows
  // Return:
  // - C: a distance matrix in a shape of (nrow(Ar), nrow(Br)).
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  mat A = mat(Ar.begin(), m, k, false);
  mat B = mat(Br.begin(), n, k, false);

  colvec An =  sum(square(A), 1);  // rowSum
  colvec Bn =  sum(square(B), 1);

  mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();

  return wrap(sqrt(C));
}

// [[Rcpp::export]]
NumericMatrix cxxcalculate_distances(NumericMatrix X){
  // Input:
  // - X: a matrix in shape of gene by sample
  // Return:
  // a distance matrix among samples
  mat X_mat = as<mat>(X);
  Rcout << "normalize to median transcript count... ";
  rowvec num_transcripts = sum(X_mat, 0);
  rowvec size_factor =  median(num_transcripts) / num_transcripts;
  NumericMatrix X_norm = wrap(X_mat.each_row() % size_factor);

  Rcout << "apply FT transformation... ";
  NumericMatrix X_FTT = cxxfreeman_tukey_transform(X_norm);

  Rcout << "calculate distance... ";
  NumericMatrix X_D = cxxdist_euclidean(transpose(X_FTT), transpose(X_FTT));

  Rcout << "[done] "<< std::endl;
  return X_D;
}
// calculate_distances <- function(mat){
//
//   num_transcripts <- Matrix::colSums(mat)
//   size_factor <- median(num_transcripts, na.rm = T) / num_transcripts
//
//   mat_norm <- t(t(mat) * size_factor)
// # apply freeman-tukey transform
//   mat_FTT <- freeman_tukey_transform(mat_norm)
// # calculate all pairwise distances using the Euclidean metric
//   mat_D <- dist(t(mat_FTT), method = "euclidean",
//                 upper = T, diag = T)
//   return(as.matrix(mat_D))
// }

/*** R
# 4 samples
x <- matrix(abs(sin(seq(from=1, to=12, length.out = 12))),
            nrow = 3, byrow = T)
print(all.equal(c(as.matrix(dist(t(x), method = "euclidean", upper=T, diag=T))),
                c(cxxdist_euclidean(t(x), t(x))) ))
print(all.equal(c(freeman_tukey_transform(x)), c(cxxfreeman_tukey_transform(x))))
print(all.equal(c(calculate_distances(x)), c(cxxcalculate_distances(x))))
*/
