#ifndef BINNED_MAT_MULT_
#define BINNED_MAT_MULT_

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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
//'
//' @export
// [[Rcpp::export]]
arma::mat binnedMatMult (const arma::mat& X, const arma::uvec& k, const arma::vec& w)
{
  unsigned int n = k.size();
  arma::mat L(X.n_cols, n, arma::fill::zeros);

  for (unsigned int i = 0; i < n; i++) {
    unsigned int ind = k(i);
    L(,ind += w(ind * X(ind,);
  }
  return L * X;
}

#endif // BINNED_MAT_MULT_
