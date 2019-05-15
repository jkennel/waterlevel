// inspired by http://davegiles.blogspot.com/2017/01/explaining-almon-distributed-lag-model.html
// and the dlnm package


#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]



arma::mat lmat(const arma::vec& x,
               int n_lags) {
  
  int n = x.n_elem;
  int n_col = n_lags;
  int n_row = n-n_col;
  
  arma::mat out = arma::mat(n_row, n_col);
  
  for (arma::uword i = 0; i < n_col; i++) {
    out.col(i) = x.subvec(i, n_row+i-1);
  }
  
  return(out);
}


// //' @title
// //' distributed_lag
// //'
// //'
// //' @description
// //' This method calculates the basis for a distributed lag.  It is currently
// //' slow.
// //'
// //' @param basisvar matrix value of lag
// //' @param basislag matrix the basis lags
// //' @param lag_max integer maximum number of lags
// //'
// //' @return distributed lag basis
// //'
// //'
// //' @export
// //'
// // [[Rcpp::export]]
// arma::mat distributed_lag_cpp(const arma::vec& bv, const arma::mat bl, arma::vec lag) {
//   
//   int nrow = bv.n_elem;
//   int ncol = bl.n_rows;
//   int max_lag = max(lag);
//   
//   arma::mat cb(ncol, nrow);
//   cb = cb.fill(NA_REAL);
//   
//   for (int v = 0; v < (nrow-max_lag); v++) {
//     cb.col(v) = arma::vectorise(bl * bv.subvec(v, v + max_lag), 0);
//   }
//   
//   return(cb.t());
//   
// }

struct dl_worker: public Worker {
  // input vector of matrices
  const arma::vec& bv;
  const arma::mat& bl;
  int lag_max;
  
  arma::mat& cb;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  dl_worker(const arma::vec& bv, const arma::mat& bl, arma::mat& cb, int lag_max)
    : bv(bv), bl(bl), cb(cb), lag_max(lag_max) {}
  
  void operator() (std::size_t begin, std::size_t end) {
    
    for (std::size_t i = begin; i < end; i++) {
      
      cb.col(i) = bl * bv.subvec(i, i + lag_max);
      
    }
  };
};





//==============================================================================
//' @title
//' distributed_lag_parallel
//'
//' @description
//' This method calculates the basis for a distributed lag in parallel.  It is currently
//' slow.
//'
//' @param x matrix value of lag
//' @param bl matrix the basis lags
//' @param lag_max integer maximum number of lags
//'
//' @return distributed lag basis
//'
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat distributed_lag_parallel(const arma::vec& x,
                                   const arma::mat& bl,
                                   int lag_max) {
  
  // result matrix
  int nrow = x.n_elem;
  int ncol = bl.n_rows;
  
  arma::mat cb(ncol, nrow);
  cb = cb.fill(NA_REAL);
  
  dl_worker calc_dl(x, bl, cb, lag_max);
  
  RcppParallel::parallelFor(0, nrow-lag_max, calc_dl);
  
  return arma::flipud(cb.t());
  
}


/***R

*/

