#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace arma;
using namespace Rcpp;


struct det_worker: public Worker {

  // input vector of matrices
  const arma::cx_cube& input;
  int n;
  arma::cx_vec& output;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  det_worker(const arma::cx_cube& input, int n, arma::cx_vec& output)
    : input(input), n(n), output(output) { }

  void operator()(std::size_t begin, std::size_t end) {

    arma::cx_mat m(n,n);

    for (std::size_t i = begin; i < end; i++) {
      m = input(arma::span(i), arma::span(), arma::span());
      output[i] = arma::det(m);

    }
  };
};


//' @title
//' det_parallel
//'
//' @description
//' Determinant for an array in parallel
//'
//' @param a \code{numeric array} values to evaluate
//'
//' @return vector of determinants
//'
//'
//' @export
// [[Rcpp::export]]
arma::cx_vec det_parallel(arma::cx_cube& a) {

  std::size_t nr = a.n_rows;
  std::size_t nc = a.n_cols;

  arma::cx_vec output(nr);

  det_worker calc_det(a, nc, output);

  RcppParallel::parallelFor(0, nr, calc_det);

  return(output);

}


//' @title
//' solve_tf_parallel
//'
//' @description
//' Calculate the transfer functions for an array in parallel
//'
//' @param a \code{numeric array} values to evaluate
//'
//' @return vector of transfer functions
//'
//'
//' @export
// [[Rcpp::export]]
arma::cx_mat solve_tf_parallel(arma::cx_cube a) {

  arma::uword nc = a.n_cols;
  arma::uword nr = a.n_rows;

  arma::cx_cube num = a(span(), span(1, nc - 1), span(1, nc - 1));
  arma::cx_vec denom = det_parallel(num);

  arma::cx_mat output(nr, nc-1);

  arma::cx_cube numer_sol = a(arma::span(),
                              arma::span(1, nc-1),
                              arma::span(0));

  for (arma::uword i = 1; i < nc; i++) {

    num(span(), span(), span(i-1)) = numer_sol;

    output.col(i - 1) = det_parallel(num) / denom;

    num(span(), span(), span(i-1)) = a(span(), span(1, nc-1), span(i));

  }


  return(output);
}

// Serial code

// //' @title
// //' det_vector
// //'
// //' @description
// //' Determinant for an array
// //'
// //' @param x \code{numeric array} values to evaluate
// //'
// //' @return vector of determinants
// //'
// //'
// //' @export
// // [[Rcpp::export]]
// arma::cx_vec det_vector(const arma::cx_cube& x) {
// 
//   std::size_t n = x.n_rows;
//   std::size_t nc = x.n_cols;
// 
//   arma::cx_vec output(n);
//   arma::cx_mat mat(nc, nc);
// 
//   for (std::size_t i = 0; i < n; i++){
//     mat = x(arma::span(i), arma::span(), arma::span());
//     output(i) = arma::det(mat);
//   }
// 
//   return(output);
// 
// }
// 
// 
// 
// // [[Rcpp::export]]
// arma::cx_mat solve_tf(arma::cx_cube x) {
// 
//   unsigned int n  = x.n_rows;
//   unsigned int nr = x.n_cols -1;
// 
//   arma::cx_mat output(n, nr);
//   arma::cx_cube numer_base = x(arma::span(),
//                                arma::span(1, nr),
//                                arma::span(1, nr));;
//   arma::cx_cube numer(n, nr, nr);
// 
// 
//   arma::cx_cube numer_sol = x(arma::span(),
//                               arma::span(1, nr),
//                               arma::span(0));
// 
//   arma::cx_vec denom = det_vector(numer_base);
// 
//   for (arma::uword i=0; i < nr; i++){
// 
//     // fill new column
//     numer = numer_base;
//     numer(arma::span(), arma::span(), arma::span(i)) = numer_sol;
// 
//     output.col(i) = det_vector(numer) / denom;
// 
//   }
// 
//   return(output);
// 
// }


/*** R

# ordering by column is the fastest
# output from spec.pgram is by row
# permuting large matrix may be slow
# currently leave output by row
# parallel achieves moderate speed-up

library(baro)
library(data.table)
library(xts)
dat <- copy(wipp30)
dat[, datetime := time * 3600 + as.POSIXct('2010-01-01')]
x <- as.matrix(dat[, list(wl, baro, et, datetime = time*3600 + as.POSIXct('1970-01-01'))])

# aa <- transfer_smooth(dat, wl = 'wl', baro = 'baro', et = 'et',
#                       datetime = 'datetime', spans = c(3,7),
#                       taper = 0.1, fast = TRUE,
#                       demean = TRUE, detrend = TRUE, pad = 0)
#

nc <- 3
nr <- 5e6
x <- matrix(rnorm(nr*nc), nrow=nr, ncol = nc)

#x <- as.matrix(dat[, list(wl, baro, et)])
tmp2 <- spec_pgram(x, spans = c(3,7),
                  taper = 0.1, fast = TRUE,
                demean = TRUE, detrend = TRUE, pad = 0, plot = FALSE)
# tmp3 <- cross_spectrum(x[,1:2], spans = c(3,7),
#                        taper = 0.1, fast = TRUE,
#                        demean = TRUE, detrend = TRUE, pad = 0)


a2 <- tmp2 * 3600
microbenchmark(
  #aa <- solve_tf(a2),
  bb <- solve_tf_parallel(a2),
  times = 2
)
all.equal(aa, bb)

*/
