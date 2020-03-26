#ifndef BINNED_MAT_MULT_
#define BINNED_MAT_MULT_

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

//' Calculating binned matrix product
//'
//' This function calculates the matrix product using Algorithm 3 of Zheyuan Li, Simon N. Wood: "Faster
//' model matrix crossproducts for large generalized linear models with discretized covariates". The idea
//' is to compute just on the unique rows of X by also using an index vector to map to the original matrix.
//' The algorithm implemented here is a small adaption of the original algorithm. Instead of calculating $XW$
//' which again, needs to be transposed, we directly calculate $X^TW$ to avoid another transposing step.
//'
//' @param X [\code{arma::mat}]\cr
//'   Matrix X.
//' @param k [\code{arma::uvec}]\cr
//'   Index vector for mapping to original matrix $X_o(i,) = X(k(i),.)$.
//' @param w [\code{arma::vec}]\cr
//'   Vector of weights that are accumulated.
//' @return \code{arma::mat} Matrix Product $X^TWX$.
//' @examples
//' nsim = 1e6L
//' nunique = trunc(sqrt(nsim))
//'
//' xunique = runif(n = nunique, min = 0, max = 10)
//' k = sample(x = seq_len(nunique), size = nsim, replace = TRUE)
//'
//' X = poly(x = xunique, degree = 20L)
//'
//' binnedMatMult(X = X, k = k-1, w = 1)
//' @export
// [[Rcpp::export]]
arma::mat binnedMatMult (const arma::mat& X, const arma::uvec& k, const arma::vec& w)
{
  unsigned int n = k.size();
  unsigned int ind;
  arma::mat L(X.n_cols, X.n_rows, arma::fill::zeros);

  // IMPROVE arma::trans(X.row(ind)): selecting X by rows is not the preferred way since
  // Armadillo uses column major layout. Maybe improve:

  if ( (w.size() == 1) && (w(0) == 1) ) {
    // std::cout << "Weight of size 1:" << std::endl;
    for (unsigned int i = 0; i < n; i++) {
       ind = k(i);
       L.col(ind) += arma::trans(X.row(ind));
    }
  } else {
    for (unsigned int i = 0; i < n; i++) {
      ind = k(i);
      L.col(ind) += arma::trans(X.row(ind)) * w(i);
    }
  }
  return L * X;
}

//' Calculating binned matrix product for resonse term
//'
//' This function calculates the matrix product using Algorithm 3 of Zheyuan Li, Simon N. Wood: "Faster
//' model matrix crossproducts for large generalized linear models with discretized covariates". The idea
//' is to compute just on the unique rows of X by also using an index vector to map to the original matrix.
//' The algorithm implemented here is a small adaption of the original algorithm. Instead of calculating $XW$
//' which again, needs to be transposed, we directly calculate $X^TW$ to avoid another transposing step. In addition
//' to the original algorithm the algorithm here directly calculates the crossproduct with the response.
//'
//' @param X [\code{arma::mat}]\cr
//'   Matrix X.
//' @param y [\code{arma::vec}]\cr
//'   Response vector y.
//' @param k [\code{arma::uvec}]\cr
//'   Index vector for mapping to original matrix $X_o(i,) = X(k(i),.)$.
//' @param w [\code{arma::vec}]\cr
//'   Vector of weights that are accumulated.
//' @return \code{arma::mat} Matrix Product $X^TWX$.
//' @examples
//' nsim = 1e6L
//' nunique = trunc(sqrt(nsim))
//'
//' xunique = runif(n = nunique, min = 0, max = 10)
//' k = sample(x = seq_len(nunique), size = nsim, replace = TRUE)
//'
//' X = poly(x = xunique, degree = 20L)
//' y = runif(nsim)
//'
//' binnedMatMultResponse(X = X, y = y, k = k-1, w = 1)
//' @export
// [[Rcpp::export]]
arma::mat binnedMatMultResponse (const arma::mat& X, const arma::vec& y,  const arma::uvec& k, const arma::vec& w)
{
  unsigned int n = k.size();
  unsigned int ind;
  arma::rowvec out(X.n_cols, arma::fill::zeros);

  // IMPROVE arma::trans(X.row(ind)): selecting X by rows is not the preferred way since
  // Armadillo uses column major layout. Maybe improve:

  if ( (w.size() == 1) && (w(0) == 1) ) {
    // std::cout << "Weight of size 1:" << std::endl;
    for (unsigned int i = 0; i < n; i++) {
       ind = k(i);
       out += y(i) * X.row(ind);
    }
  } else {
    for (unsigned int i = 0; i < n; i++) {
      ind = k(i);
      out += w(i) * y(i) * X.row(ind);
    }
  }
  return out;
}

#endif // BINNED_MAT_MULT_
