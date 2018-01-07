#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]
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
NumericMatrix cxxcalculate_distances(NumericMatrix X, int verbose=0){
  // Input:
  // - X: a matrix in shape of gene by sample
  // - verbose: verbose level
  // Return:
  // a distance matrix among samples
  mat X_mat = as<mat>(X);
  if (verbose > 3) Rcout << "normalize to median transcript count... ";
  rowvec num_transcripts = sum(X_mat, 0);
  rowvec size_factor =  median(num_transcripts) / num_transcripts;
  NumericMatrix X_norm = wrap(X_mat.each_row() % size_factor);

  if (verbose > 3) Rcout << "apply FT transformation... ";
  NumericMatrix X_FTT = cxxfreeman_tukey_transform(X_norm);

  if (verbose > 3) Rcout << "calculate distance... ";
  NumericMatrix X_D = cxxdist_euclidean(transpose(X_FTT), transpose(X_FTT));

  if (verbose > 3) Rcout << "[done] "<< std::endl;
  return X_D;
}

// [[Rcpp::export]]
NumericMatrix cxxknn_smoothing_onestep(NumericMatrix Dr,
                                       NumericMatrix Xr,
                                       int k=2,
                                       int verbose=0){
  // Input:
  // - D: a distance matrix in a shape of #samples x #samples
  // - X: a matrix in a shape of #genes x #samples
  // - k: an integer specifying the k neighbouring samples (self-inclusive) to choose
  // Return:
  // - S: a matrix in a same shape of X after summing up k neighbours.
  if (k > Xr.ncol()) stop("k is greater than the number of samples.");

  if (verbose > 3) Rcout << "aggregating neighbours... " ;
  mat D = as<mat>(Dr);
  mat X = as<mat>(Xr);
  mat S(X.n_rows, X.n_cols);

  // umat A = D.each_col( [](vec& a){ sort_index(a); } ); // lambda-func fail
  umat A(k, D.n_cols);
  for (int j = 0; j < A.n_cols; j++){
    uvec idx = sort_index(D.col(j));
    uvec idx_k = idx.head(k);
    A.col(j) = idx_k;
  }
  for (int j = 0; j < S.n_cols; j++){
    uvec idx_k_nb = A.col(j);
    S.col(j) = sum(X.cols(idx_k_nb), 1);
  }
  if (verbose > 3) Rcout << "[done]" << std::endl;
  return wrap(S);
}

// [[Rcpp::export]]
NumericMatrix cxxknn_smoothing(NumericMatrix Xr,
                               int k=5,
                               int verbose=0){
  if (k > Xr.ncol()) stop("k is greater than the number of samples.");
  NumericMatrix S = clone(Xr);
  int num_powers = ceil(log1p(k) / log(2));

  for (int p = 1; p < num_powers + 1; p++){

    int k_step = int(pow(2, p)-1);
    if (k_step > k) k_step = k;

    if (verbose > 3) Rcout << "Step " << p << "/" << num_powers << ": " << \
        "Smoothing with k=" << k_step << " " << std::endl;
    NumericMatrix D = cxxcalculate_distances(S, verbose);
    S = cxxknn_smoothing_onestep(D, Xr, k_step+1, verbose);
  }
  if (verbose > 3) Rcout << "[KNN-Smoothing Finished]" << std::endl;
  return S;
}
/*** R
num_genes = 5
num_samples = 10
x <- matrix(abs(sin(seq(from=1, to=num_samples*num_genes))) * 100,
            ncol = num_samples, byrow = T)
colnames(x) <- paste0('cell', seq_len(ncol(x)))
rownames(x) <- paste0('gene', seq_len(nrow(x)))
xd <- as.matrix(dist(t(x)))
print(all.equal(c(as.matrix(dist(t(x), method = "euclidean", upper=T, diag=T))),
                c(cxxdist_euclidean(t(x), t(x))) ))
print(all.equal(c(freeman_tukey_transform(x)), c(cxxfreeman_tukey_transform(x))))
print(all.equal(c(calculate_distances(x)), c(cxxcalculate_distances(x))))
print(all.equal(c(knn_smoothing(x)), c(cxxknn_smoothing(x, verbose=4))))
*/
