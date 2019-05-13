// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

//' @title
//' shift_subset
//'
//' @description
//' lag data and subset the results
//'
//' @param x \code{numeric vector} to lag
//' @param lag \code{integer vector} with the lags
//' @param n_subset \code{integer} subset every n_subset values
//' @param n_shift \code{integer} amount to shift results
//'
//' @return vector with lagged values
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector shift_subset(Rcpp::NumericVector x,
                                 int lag,
                                 int n_subset,
                                 int n_shift) {
  
  int n = x.size();
  int first = n_shift + lag;
  int n_values;
  int k;
  int wh;
  int end_adj;
  
  if(n_subset == 1){
    n_values = (n-n_shift);
  } else {
    n_values = ((n-n_shift-1) / n_subset) + 1;
  }
  
  if (lag > n_shift) {
    k = 1;
    k = k + (lag-1) / (n_subset);
  } else {
    k = 0;
    end_adj = ((-lag + n_shift) / n_subset);
    if ((-lag + n_shift) % n_subset != 0) {
      end_adj = end_adj + 1;
    }
    if (((-lag + n_shift) / n_subset + 1))
      n_values = n_values - end_adj;
  }
  
  Rcpp::NumericVector out(n, NA_REAL);
  
  for (std::size_t i = k; i < n_values; i++) {
    
    wh = i * n_subset + n_shift - lag;
    
    if (wh < n) {
      
      out[i] = x[wh];
      
    } else {
      throw std::range_error("lag is too large");
      // Rcout << "lag is " << lag << std::endl;
      // Rcout << "wh is " << wh << std::endl;
    }
  }
  //}
  
  return out;
}



//' @title
//' lag_matrix
//'
//' @description
//' lag data and subset the results
//'
//' @param x \code{numeric vector} to lag
//' @param lags \code{integer vector} with the lags
//' @param n_subset \code{integer} subset every n_subset values
//' @param n_shift \code{integer} amount to shift results
//' @param var_name \code{character} name for the generated matrix columns
//'
//' @return matrix with lagged values
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix lag_matrix(Rcpp::NumericVector x,
                               Rcpp::IntegerVector lags,
                               int n_subset = 1,
                               int n_shift = 0,
                               std::string var_name = "lag") {
  
  int n = x.size();
  int n_row;
  
  if(n_subset == 1){
    n_row = (n-n_shift);
  } else {
    n_row = ((n-n_shift-1) / n_subset) + 1;
  }
  int n_cols = lags.size();
  
  Rcpp::CharacterVector nm(n_cols);
  Rcpp::NumericMatrix out = Rcpp::NumericMatrix(n_row, n_cols);
  
  for (std::size_t i = 0; i < n_cols; i++) {
    out(_, i) = shift_subset(x, lags[i], n_subset, n_shift);
    if(lags[i] < 0) {
      nm[i] = var_name + '_' + 'n' + std::to_string(abs(lags[i]));
    } else {
      nm[i] = var_name + '_' + std::to_string(lags[i]);
    }
  }
  
  colnames(out) = nm;
  
  return(out);
}
