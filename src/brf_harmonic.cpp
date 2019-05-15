// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// Armadillo does sin cos much faster than base R.  Together there is also some 
// optimization.

//==============================================================================
//' harmonic
//'
//' calculate sin and cos curves from POSIXct times (serial)
//'
//' @param times \code{numeric vector} times to calculate sin and cos at
//' @param freq \code{numeric vector} frequencies for sin and cos
//'
//' @return sin and cos curves
//'
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix harmonic(arma::vec times,
                             arma::vec& freq){
  
  unsigned int n_rows = times.n_elem;
  unsigned int n_cols = freq.n_elem;
  times -= times(0);
  Rcpp::CharacterVector nm(n_cols * 2);
  
  
  arma::mat out = arma::mat(n_rows, n_cols * 2);
  
  // 2 * pi * times / 86400
  arma::vec cycle = (M_PI * (times / 43200));
  
  for (std::size_t i = 0; i < n_cols; i++) {
    
    //harmonic = cycle * freq(i);
    out.col(i)          = arma::sin(cycle * freq(i));
    out.col(i + n_cols) = arma::cos(cycle * freq(i));
    
  }
  
  for (std::size_t i = 0; i < n_cols; i++) {
    nm[i]          = "sin_" + std::to_string(freq[i]);
    nm[i + n_cols] = "cos_" + std::to_string(freq[i]);
  }
  
  
  // add names
  Rcpp::NumericMatrix out_(wrap(out)); 
  Rcpp::colnames(out_) = nm;
  
  
  return(out_);
}




// //==============================================================================
// //' sin_harmonic
// //'
// //' calculate sin and cos curves from POSIXct times (serial)
// //'
// //' @param x \code{numeric vector} times to calculate sin and cos at
// //' @param freq \code{numeric vector} frequencies for sin and cos
// //'
// //' @return sin curve
// //'
// //'
// //' @export
// //'
// // [[Rcpp::export]]
// arma::vec sin_harmonic(arma::vec x,
//                                  int freq){
//   
//   unsigned int n_rows = x.n_elem;
//   x -= x(0);
// 
//   arma::vec out = arma::vec(n_rows);
//   
//   // 2 * pi * x / 86400
//   arma::vec cycle = (M_PI * (x / 43200));
//   
// 
//   out = arma::sin(cycle * freq);
//   
//   return(out);
//   
// }
// 
// 
// //==============================================================================
// //' cos_harmonic
// //'
// //' calculate sin and cos curves from POSIXct times (serial)
// //'
// //' @param x \code{numeric vector} times to calculate sin and cos at
// //' @param freq \code{numeric vector} frequencies for sin and cos
// //'
// //' @return cos curve
// //'
// //' @export
// //'
// // [[Rcpp::export]]
// arma::vec cos_harmonic(arma::vec x,
//                        int freq){
//   
//   int n_rows = x.n_elem;
//   x -= x(0);
//   
//   arma::vec out = arma::vec(n_rows);
//   
//   // 2 * pi * x / 86400
//   arma::vec cycle = (M_PI * (x / 43200));
//   
//   
//   out = arma::cos(cycle * freq);
//   
//   return(out);
//   
// }

