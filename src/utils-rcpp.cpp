#include <RcppArmadillo.h>
#include <bigmemory/MatrixAccessor.hpp>
// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



// Logic for BigRowSums.
template <typename T>
arma::vec BigRowSums(XPtr<BigMatrix> pMat, MatrixAccessor<T> mat) {
  arma::vec rowSums(pMat->nrow(), fill::zeros);
  arma::vec value(1);
  for (int jj = 0; jj < pMat->ncol(); jj++) {
    for (int ii = 0; ii < pMat->nrow(); ii++) {
      value = mat[jj][ii];
      if (!value.has_nan()) {
        rowSums[ii] += value[0];
      }
    }
  }
  return rowSums;
}
// Dispatch function for BigRowSums
//
// [[Rcpp::export]]
arma::vec BigRowSums(SEXP pBigMat) {
  XPtr<BigMatrix> xpMat(pBigMat);

  switch(xpMat->matrix_type()) {
  case 1:
    return BigRowSums(xpMat, MatrixAccessor<char>(*xpMat));
  case 2:
    return BigRowSums(xpMat, MatrixAccessor<short>(*xpMat));
  case 4:
    return BigRowSums(xpMat, MatrixAccessor<int>(*xpMat));
  case 6:
    return BigRowSums(xpMat, MatrixAccessor<float>(*xpMat));
  case 8:
    return BigRowSums(xpMat, MatrixAccessor<double>(*xpMat));
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}


// 
// 
// // [[Rcpp::export]]
// vec uni_WQR_by_col(mat& Y, vec& x, vec& w, const double tau, const bool rm_bas) {
//   
//   // DESCRIPTION
//   // Function to compute the univariate weighted quantile regression over the column of 
//   // a matrix Y given the covariate x, the weights w, and the percentile to fit tau.
//   
//   // INPUT
//   // Y: matrix where each column will be regressed on
//   // x: the covariate
//   // w: vector of weights associated to rows of Y 
//   // tau: the percentile to fit 
//   
//   // OUTPUT
//   // b: a vector containing the fitted coefs for every column of Y
//   
//   
//   // Subset the indices for which the covariate and the weights are not zero
//   // uvec nz = find(x  > 0);
//   uvec nz = find(x % w > 0);
//   vec x_nz = x.elem(nz);
//   vec w_nz = w.elem(nz);
//   mat Y_nz = Y.rows(nz);
//   // Compute the weight threshold that will be used to select the index
//   const double thres = tau * sum(w_nz);
//   // Initialize the vector of coefficients beta
//   const int c = Y_nz.n_cols;
//   vec b(c);
//   int i;
//   // Iterate over the columns of Y
//   for (i = 0; i < c; i++) {
//     // Compute the vector of ratios y/x
//     vec y_nz = Y_nz.col(i); 
//     if(rm_bas){
//       y_nz = y_nz - min(y_nz);
//     } 
//     vec r = y_nz / x_nz;
//     // Sort the vector of ratios and the vector of weights
//     uvec idx = sort_index(r, "ascend");
//     vec r_sorted = r(idx);
//     vec w_sorted = w_nz(idx);
//     // Compute the criterion vector
//     vec w_crit = cumsum(w_sorted) - thres;
//     // Find the index where the criterion vector becomes non-negative
//     unsigned int j; // unsigned because j will never become negative
//     for(j = 0; j < w_crit.size(); j++){
//       if(w_crit(j) >= 0){ 
//         break;
//       }
//     }
//     // Store coefficient in output vector
//     b(i) = r_sorted(j-1);
//   }
//   return b;
// }