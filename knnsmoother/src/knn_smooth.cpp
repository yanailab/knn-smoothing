#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]
using namespace Rcpp;
using namespace arma;


//' Freeman-Tukey transformation on a matrix
//'
//' \code{sqrt(X) + sqrt(X + 1)}
//'
//' @param X A matrix in a shape of #genes x #samples.
//' @return A matrix in a same shape of \code{X} with freeman tukey
//' transformation performed.
//' @export
// [[Rcpp::export]]
NumericMatrix freeman_tukey_transform(NumericMatrix X){
  mat X_mat = as<mat>(X);
  return wrap(sqrt(X_mat) + sqrt(X_mat + 1));
}


//' Pair-wise euclidean distance among two matrices.
//'
//' @param Ar A matrix with samples on rows.
//' @param Br A matrix with samples on rows.
//' @return A distance matrix in a shape of nrow(\code{Ar}) x nrow(\code{Br}).
//' @references \url{http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/}
//' @examples
//' x <- matrix(1:12, 3)
//' all.equal(c(as.matrix(dist(x))), c(dist_euclidean(x, x)))
//' @export
// [[Rcpp::export]]
NumericMatrix dist_euclidean(NumericMatrix Ar, NumericMatrix Br) {

  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  mat A = mat(Ar.begin(), m, k, false);
  mat B = mat(Br.begin(), n, k, false);

  colvec An =  sum(square(A), 1);  // rowSum
  colvec Bn =  sum(square(B), 1);

  mat C = -2 * (A * B.t());
  C.each_col() += An; checkUserInterrupt();
  C.each_row() += Bn.t(); checkUserInterrupt();

  return wrap(sqrt(C));
}


//' Distance matrix
//'
//' Normalize (by median) and apply Freeman-Tukey transformation on the input
//' matrix \code{X}. Then calculate the distance matrix of samples.
//'
//' @param X A matrix in a shape of #genes x #samples.
//' @param verbose An integer to specify verbose level.
//' @return A distance matrix in a shape of #samples x #samples.
//' @export
// [[Rcpp::export]]
NumericMatrix smoother_calc_distance(NumericMatrix X, int verbose=0){
  mat X_mat = as<mat>(X);
  if (verbose > 2) Rcout << "normalize to median transcript count... ";
  rowvec num_transcripts = sum(X_mat, 0);
  rowvec size_factor =  median(num_transcripts) / num_transcripts;
  NumericMatrix X_norm = wrap(X_mat.each_row() % size_factor);

  if (verbose > 2) Rcout << "apply FT transformation... ";
  NumericMatrix X_FTT = freeman_tukey_transform(X_norm);

  if (verbose > 2) Rcout << "calculate distance... ";
  NumericMatrix X_D = dist_euclidean(transpose(X_FTT), transpose(X_FTT));

  if (verbose > 2) Rcout << "[done] "<< std::endl;
  return X_D;
}


//' Aggregate K nearest expression profiles
//'
//' Normalize (by median) and apply Freeman-Tukey transformation on the input
//' matrix \code{X}. Then calculate the distance matrix of samples.
//'
//' @param Xr A matrix in a shape of #genes x #samples.
//' @param Dr A predefined distance matrix in a shape of #samples x #samples. If
//'   not specified, D is the distance matrix of the input \code{Xr}.
//' @param k An integer to choose \code{k} nearest samples (self-inclusive) to
//'  aggregate based on the distance matrix \code{Dr}. If \code{k} is greater than
//'  #samples, \code{k} is forced to be #samples to continue aggregation.
//' @param verbose An integer to specify verbose level.
//' @return An aggregated matrix in a same shape of \code{Xr}.
//' @export
// [[Rcpp::export]]
NumericMatrix aggregate_k_nearest(
    NumericMatrix Xr,
    NumericMatrix Dr = NumericMatrix(),
    int k=2,
    int verbose=0) {
  if (k > Xr.ncol()) {
    Rcout << "Specified k=" << k << \
      " is greater than #samples=" << Xr.ncol() << ". " << \
      "Force k=#samples." << std::endl;
    k = Xr.ncol();
  }
  if (all(is_na(Dr))) Dr = smoother_calc_distance(Xr);

  if (verbose > 2) Rcout << "aggregating neighbours... " ;
  mat D = as<mat>(Dr);
  mat X = as<mat>(Xr);
  mat S(X.n_rows, X.n_cols);

  // umat A = D.each_col( [](vec& a){ sort_index(a); } ); // lambda-func fails
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
  if (verbose > 2) Rcout << "[done]" << std::endl;
  return wrap(S);
}


//' Perform KNN-smoothing on UMI-filtered scRNA-seq data
//'
//' @param X A matrix in a shape of #genes x #samples.
//' @param k An integer to choose \code{k} nearest samples (self-inclusive) to
//'  aggregate based on the distance matrix \code{Dr}. If \code{k} is greater than
//'  #samples, \code{k} is forced to be #samples to continue aggregation.
//' @param verbose An integer to specify verbose level.
//' @return An aggregated matrix in a same shape of \code{X}.
//' @references  "K-nearest neighbor smoothing for high-throughput single-cell
//'   RNA-Seq data" (Florian Wagner, Yun Yan, Itai Yanai, bioRxiv 217737; doi:
//'   \url{https://doi.org/10.1101/217737}).
//' @export
// [[Rcpp::export]]
NumericMatrix knn_smoothing(NumericMatrix X,
                            int k=5,
                            int verbose=0){
  if (k > X.ncol()) stop("k should not be greater than #samples.");
  NumericMatrix S = clone(X);

  int num_powers = ceil(log1p(k) / log(2));
  for (int p = 1; p < num_powers + 1; p++){
    int k_step = int(pow(2, p)-1);
    if (k_step > k) k_step = k;

    if (verbose > 2) Rcout << "Step " << p << "/" << num_powers << ": " << \
        "Smoothing with k=" << k_step << " " << std::endl;

    NumericMatrix D = smoother_calc_distance(S, verbose);
    checkUserInterrupt();
    S = aggregate_k_nearest(X, D, k_step+1, verbose);
    checkUserInterrupt();
  }
  if (verbose > 2) Rcout << "[KNN-Smoothing Finished]" << std::endl;
  return S;
}

