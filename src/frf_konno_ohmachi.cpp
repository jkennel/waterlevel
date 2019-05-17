//Adapted from https://github.com/jsh9/fast-konno-ohmachi/blob/master/konno_ohmachi.py

// Copyright (c) 2018 Jonathan Kennel
// MIT License
// 
// Copyright (c) 2013-2017, Jian Shi
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//   
//   The above copyright notice and this permission notice shall be included in all
//   copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//   SOFTWARE.

#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// x spectra vals
// f frequency
// b smoothing coef (even)

// [[Rcpp::export]]
arma::vec calc_b_vals(int b,
                      const arma::vec& wb) {

  arma::vec ret = arma::pow((sin(b * log10(wb)) / (b * log10(wb))), 4);
  ret.elem(arma::find_nonfinite(ret)).ones();

  return(ret);
}

// [[Rcpp::export]]
double konno_ohmachi(const arma::vec& b_vals,
                     const arma::vec& ref_z,
                     const arma::vec& f,
                     const arma::vec& x,
                     int i) {

  double fc, sum_w;
  int n = x.n_elem;
  int n_w0;
  int idx, s, e;

  arma::uvec ind(n);
  arma::vec w(n);
  w.zeros();

  fc  = f[i];
  ind = arma::find(f >= 0.5 * fc && f <= 2.0 * fc);

  arma::vec z = f(ind) / fc;

  n_w0 = z.n_elem;

  arma::vec w0(n_w0);
  arma::interp1(ref_z, b_vals, z, w0, "*linear", 1.0);

  // arma::vec w0 = arma::pow((sin(10 * log10(z)) / (10 * log10(z))), 4);
  // w0.elem(arma::find_nonfinite(w0)).zeros();

  idx = arma::index_max(w0);
  s = (i - idx);

  // shift are values to pad
  if ((s + n_w0) > n) {
    e = n - 1;
    w.subvec(s, e) = w0;
    sum_w = arma::sum(w);
  } else {
    e = s + n_w0 - 1;
    w.subvec(s, e) = w0;
    sum_w = arma::sum(w0);
  }

  return(arma::dot(w, x) / sum_w);

}


struct ko_worker: public Worker {
  // input vector of matrices
  const arma::vec& b_vals;
  const arma::vec& ref_z;
  const arma::vec& f;
  const arma::vec& x;
  arma::vec& y;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  ko_worker(const arma::vec& b_vals,
            const arma::vec& ref_z,
            const arma::vec& f,
            const arma::vec& x,
            arma::vec& y)
    : b_vals(b_vals), ref_z(ref_z), f(f), x(x), y(y) { }

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t i = begin; i < end; i++) {

      y(i) = konno_ohmachi(b_vals, ref_z, f, x, i);

    }
  };
};


//==============================================================================
//' @title
//' konno_ohmachi_parallel
//'
//' @description
//' This method does konno ohmachi smoothing
//'
//' @param x vector to smooth
//' @param f vector of frequencies
//' @param b integer even magnitude to smooth
//'
//' @return konno ohmachi smoothed results
//'
//'
//' @export
//'
// [[Rcpp::export]]
arma::vec konno_ohmachi_parallel(const arma::vec& x,
                                 const arma::vec& f,
                                 int b = 10) {

  int n = x.n_elem;
  arma::vec y(n);
  y.zeros();

  arma::vec ref_z = arma::regspace<arma::vec>(0.5, 0.0005, 2.0);
  arma::vec b_vals = calc_b_vals(b, ref_z);

  ko_worker calc_ko(b_vals, ref_z, f, x, y);


  RcppParallel::parallelFor(1, n, calc_ko);


  y(0) = y(1);
  y(n-1) = y(n-2);


  return(y);
}


//==============================================================================
//' @title
//' konno_ohmachi_serial
//'
//' @description
//' This method does konno ohmachi smoothing
//'
//' @param x \code{numeric vector} to smooth
//' @param f \code{numeric vector} of frequencies
//' @param b \code{integer} smoothing coefficient
//'
//' @return vector of konno ohmachi smoothed results
//'
//'
//' @export
//'
// [[Rcpp::export]]
arma::vec konno_ohmachi_serial(const arma::vec& x,
                               const arma::vec& f,
                               int b = 10) {

  int n = x.n_elem;
  arma::vec y(n);
  y.zeros();

  arma::vec ref_z = arma::regspace<arma::vec>(0.5, 0.0005, 2.0);
  arma::vec b_vals = calc_b_vals(b, ref_z);

  for (int i = 1; i < n; i++) {
    
    y(i) = konno_ohmachi(b_vals, ref_z, f, x, i);

  }

  y(0) = y(1);
  y(n-1) = y(n-2);


  return(y);
}




